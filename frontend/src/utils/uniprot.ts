const ACCESSION_RE = /^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\d+)?$/i

/**
 * Extract a UniProt accession from various formats:
 *   sp|A5FN12|ENO_FLAJ1  -> A5FN12
 *   tr|A0A4X1U3H6|A0A4X1U3H6_PIG -> A0A4X1U3H6
 *   A5FN12               -> A5FN12
 *   A5FN12-2             -> A5FN12-2 (isoform)
 * Returns null if no valid accession is found.
 */
export function extractUniprotAccession(id: string): string | null {
  const s = id.trim()

  // Pipe-delimited: db|accession|entry_name
  const parts = s.split('|')
  if (parts.length >= 2) {
    const db = parts[0].toLowerCase()
    if (db === 'sp' || db === 'tr') {
      const acc = parts[1]
      if (ACCESSION_RE.test(acc)) return acc
    }
  }

  // Bare accession (with optional isoform suffix)
  if (ACCESSION_RE.test(s)) return s

  return null
}

export function isUniprotAccession(id: string): boolean {
  return extractUniprotAccession(id) !== null
}

export function uniprotUrl(accession: string): string {
  const base = accession.split('-')[0].toUpperCase()
  return `https://www.uniprot.org/uniprotkb/${base}`
}

export interface UniprotInfo {
  accession: string
  name: string
  fullName: string
  organism: string
  gene: string | null
  reviewed: boolean
}

const cache = new Map<string, UniprotInfo | null>()
const inFlight = new Map<string, Promise<UniprotInfo | null>>()

export async function fetchUniprotInfo(accession: string): Promise<UniprotInfo | null> {
  const key = accession.split('-')[0].toUpperCase()

  if (cache.has(key)) return cache.get(key)!

  if (inFlight.has(key)) return inFlight.get(key)!

  const promise = (async (): Promise<UniprotInfo | null> => {
    try {
      const res = await fetch(
        `https://rest.uniprot.org/uniprotkb/${key}?format=json`,
        { signal: AbortSignal.timeout(8000) }
      )
      if (!res.ok) {
        cache.set(key, null)
        return null
      }
      const data = await res.json()

      const recommendedName =
        data.proteinDescription?.recommendedName?.fullName?.value ??
        data.proteinDescription?.submissionNames?.[0]?.fullName?.value ??
        null

      const shortName =
        data.proteinDescription?.recommendedName?.shortNames?.[0]?.value ?? null

      const organism: string =
        data.organism?.scientificName ?? data.organism?.commonName ?? ''

      const gene: string | null =
        data.genes?.[0]?.geneName?.value ?? null

      const reviewed: boolean = data.entryType === 'UniProtKB reviewed (Swiss-Prot)'

      const info: UniprotInfo = {
        accession: key,
        name: shortName ?? recommendedName ?? key,
        fullName: recommendedName ?? key,
        organism,
        gene,
        reviewed,
      }

      cache.set(key, info)
      return info
    } catch {
      cache.set(key, null)
      return null
    } finally {
      inFlight.delete(key)
    }
  })()

  inFlight.set(key, promise)
  return promise
}
