import { useEffect, useState } from 'react'
import { Loader2, ChevronRight, ChevronDown, Download } from 'lucide-react'
import { registerMappingFile } from '../utils/duckdb'
import { getDuckDB } from '../utils/duckdb'
import { isUniprotAccession } from '../utils/uniprot'
import UniprotProteinLabel from './UniprotProteinLabel'

interface MappingRow {
  peptide: string
  peptide_lca_tax_ids: number[]
  peptide_go_terms: string[]
  background_protein: string
  annotated_protein: string
  evalue: number | null
  pident: number | null
}

interface SubjectHit {
  evalue: number | null
  pident: number | null
}

type BgProteinMap = Map<string, Map<string, SubjectHit>>
type PeptideMap = Map<string, BgProteinMap>

interface PeptideDetailsPaneProps {
  jobId: string
  listId: string
  selectedTaxIds?: string[] | null
  selectedTaxName?: string | null
  selectedGoId?: string | null
  selectedGoName?: string | null
}

export default function PeptideDetailsPane({
  jobId,
  listId,
  selectedTaxIds,
  selectedTaxName,
  selectedGoId,
  selectedGoName,
}: PeptideDetailsPaneProps) {
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState<string | null>(null)
  const [hierarchy, setHierarchy] = useState<PeptideMap>(new Map())
  const [expandedPeptides, setExpandedPeptides] = useState<Set<string>>(new Set())
  const [expandedBgProteins, setExpandedBgProteins] = useState<Set<string>>(new Set())
  // registrationKey tracks which job/list the view is registered for
  const [registrationKey, setRegistrationKey] = useState<string | null>(null)

  const currentKey = `${jobId}/${listId}`

  useEffect(() => {
    if (registrationKey === currentKey) return
    setRegistrationKey(null)

    registerMappingFile(jobId, listId)
      .then(() => {
        setRegistrationKey(currentKey)
      })
      .catch((err) => {
        setError(err instanceof Error ? err.message : 'Failed to initialize data engine')
      })
  }, [jobId, listId, currentKey, registrationKey])

  useEffect(() => {
    // Don't query until the view for this job/list is registered
    if (registrationKey !== currentKey) return

    const hasTax = selectedTaxIds && selectedTaxIds.length > 0
    if (!hasTax && !selectedGoId) {
      setHierarchy(new Map())
      setExpandedPeptides(new Set())
      setExpandedBgProteins(new Set())
      return
    }

    let cancelled = false

    async function runQuery() {
      setLoading(true)
      setError(null)
      try {
        const { conn } = await getDuckDB()

        const conditions: string[] = []
        if (selectedTaxIds && selectedTaxIds.length > 0) {
          const ids = selectedTaxIds.map(id => parseInt(id, 10)).filter(n => !isNaN(n))
          if (ids.length > 0) {
            conditions.push(`list_has_any(peptide_lca_tax_ids, [${ids.join(', ')}])`)
          }
        }
        if (selectedGoId) {
          const escaped = selectedGoId.replace(/'/g, "''")
          conditions.push(`list_contains(peptide_go_terms, '${escaped}')`)
        }

        const where = conditions.length > 0 ? `WHERE ${conditions.join(' AND ')}` : ''
        const sql = `
          SELECT peptide, peptide_lca_tax_ids, peptide_go_terms, background_protein, annotated_protein, evalue, pident
          FROM mappings
          ${where}
        `

        const result = await conn.query(sql)
        if (cancelled) return

        const map: PeptideMap = new Map()
        const rows = result.toArray()
        for (const row of rows) {
          const r = row.toJSON() as MappingRow
          if (!map.has(r.peptide)) map.set(r.peptide, new Map())
          const bgMap = map.get(r.peptide)!
          if (!bgMap.has(r.background_protein)) bgMap.set(r.background_protein, new Map())
          bgMap.get(r.background_protein)!.set(r.annotated_protein, { evalue: r.evalue ?? null, pident: r.pident ?? null })
        }

        setHierarchy(map)
        setExpandedPeptides(new Set())
        setExpandedBgProteins(new Set())
      } catch (err) {
        if (!cancelled) {
          setError(err instanceof Error ? err.message : 'Query failed')
        }
      } finally {
        if (!cancelled) setLoading(false)
      }
    }

    runQuery()
    return () => { cancelled = true }
  }, [selectedTaxIds, selectedGoId, registrationKey, currentKey])

  function togglePeptide(peptide: string) {
    setExpandedPeptides(prev => {
      const next = new Set(prev)
      if (next.has(peptide)) next.delete(peptide)
      else next.add(peptide)
      return next
    })
  }

  function toggleBgProtein(key: string) {
    setExpandedBgProteins(prev => {
      const next = new Set(prev)
      if (next.has(key)) next.delete(key)
      else next.add(key)
      return next
    })
  }

  function downloadCsv() {
    const escape = (s: string) => `"${s.replace(/"/g, '""')}"`
    const lines: string[] = ['peptide,background_protein,annotated_protein,evalue,pident']
    for (const [peptide, bgMap] of Array.from(hierarchy.entries()).sort(([a], [b]) => a.localeCompare(b))) {
      for (const [bgProtein, subjects] of Array.from(bgMap.entries())) {
        for (const [subject, hit] of Array.from(subjects.entries()).sort(([a], [b]) => a.localeCompare(b))) {
          const evalue = hit.evalue !== null ? String(hit.evalue) : ''
          const pident = hit.pident !== null ? String(hit.pident) : ''
          lines.push(`${escape(peptide)},${escape(bgProtein)},${escape(subject)},${evalue},${pident}`)
        }
      }
    }
    const blob = new Blob([lines.join('\n')], { type: 'text/csv' })
    const url = URL.createObjectURL(blob)
    const a = document.createElement('a')
    a.href = url
    a.download = 'peptide_details.csv'
    a.click()
    URL.revokeObjectURL(url)
  }

  const hasSelection = (selectedTaxIds && selectedTaxIds.length > 0) || selectedGoId

  return (
    <div className="flex flex-col h-full border border-gray-200 dark:border-gray-700 rounded-lg bg-white dark:bg-gray-900 overflow-hidden">
      <div className="px-4 py-3 border-b border-gray-100 dark:border-gray-700 bg-gray-50 dark:bg-gray-800">
        <div className="flex items-center justify-between">
          <h3 className="text-sm font-semibold text-gray-700 dark:text-gray-200">Peptide Details</h3>
          {hasSelection && !loading && hierarchy.size > 0 && (
            <button
              onClick={downloadCsv}
              title="Download as CSV"
              className="p-1 rounded text-gray-400 dark:text-gray-500 hover:text-gray-700 dark:hover:text-gray-200 hover:bg-gray-200 dark:hover:bg-gray-700 transition-colors"
            >
              <Download className="w-3.5 h-3.5" />
            </button>
          )}
        </div>
        {hasSelection && (
          <p className="text-xs text-gray-500 dark:text-gray-400 mt-0.5">
            {selectedTaxIds && selectedTaxIds.length > 0 && (
              <span>{selectedTaxName ?? `Tax ID: ${selectedTaxIds[0]}`}</span>
            )}
            {selectedTaxIds && selectedTaxIds.length > 0 && selectedGoId && <span> · </span>}
            {selectedGoId && (
              <span>{selectedGoName ? `${selectedGoId} (${selectedGoName})` : selectedGoId}</span>
            )}
          </p>
        )}
      </div>

      <div className="flex-1 overflow-y-auto">
        {!hasSelection && (
          <div className="flex items-center justify-center h-full text-gray-400 dark:text-gray-500 text-sm px-4 text-center">
            Click a node in the visualization to see peptide details.
          </div>
        )}

        {hasSelection && loading && (
          <div className="flex items-center justify-center h-full gap-2 text-gray-500 dark:text-gray-400 text-sm">
            <Loader2 className="w-4 h-4 animate-spin" />
            Loading…
          </div>
        )}

        {hasSelection && !loading && error && (
          <div className="p-4 text-sm text-red-600 dark:text-red-400">{error}</div>
        )}

        {hasSelection && !loading && !error && hierarchy.size === 0 && (
          <div className="flex items-center justify-center h-full text-gray-400 dark:text-gray-500 text-sm px-4 text-center">
            No peptides found for this selection.
          </div>
        )}

        {hasSelection && !loading && !error && hierarchy.size > 0 && (
          <ul className="divide-y divide-gray-100 dark:divide-gray-800">
            {Array.from(hierarchy.entries()).sort(([a], [b]) => a.localeCompare(b)).map(([peptide, bgMap]) => {
              const isExpanded = expandedPeptides.has(peptide)
              return (
                <li key={peptide}>
                  <button
                    className="w-full flex items-center gap-2 px-4 py-2 text-left hover:bg-gray-50 dark:hover:bg-gray-800 transition-colors"
                    onClick={() => togglePeptide(peptide)}
                  >
                    {isExpanded
                      ? <ChevronDown className="w-3.5 h-3.5 text-gray-400 flex-shrink-0" />
                      : <ChevronRight className="w-3.5 h-3.5 text-gray-400 flex-shrink-0" />}
                    <span className="font-mono text-xs text-gray-800 dark:text-gray-200 break-all">{peptide}</span>
                    <span className="ml-auto text-xs text-gray-400 dark:text-gray-500 flex-shrink-0">{bgMap.size} prot.</span>
                  </button>

                  {isExpanded && (
                    <ul className="bg-gray-50 dark:bg-gray-800/50 border-t border-gray-100 dark:border-gray-700">
                      {Array.from(bgMap.entries()).map(([bgProtein, subjects]) => {
                        const bgKey = `${peptide}::${bgProtein}`
                        const isBgExpanded = expandedBgProteins.has(bgKey)
                        return (
                          <li key={bgProtein}>
                            <button
                              className="w-full flex items-center gap-2 pl-8 pr-4 py-1.5 text-left hover:bg-gray-100 dark:hover:bg-gray-700 transition-colors"
                              onClick={() => toggleBgProtein(bgKey)}
                            >
                              {isBgExpanded
                                ? <ChevronDown className="w-3 h-3 text-gray-400 flex-shrink-0" />
                                : <ChevronRight className="w-3 h-3 text-gray-400 flex-shrink-0" />}
                              <span className="font-mono text-xs text-indigo-700 dark:text-indigo-400 break-all">{bgProtein}</span>
                              <span className="ml-auto text-xs text-gray-400 dark:text-gray-500 flex-shrink-0">{subjects.size} subj.</span>
                            </button>

                            {isBgExpanded && (
                              <ul className="bg-white dark:bg-gray-900 border-t border-gray-100 dark:border-gray-700">
                                {Array.from(subjects.entries()).sort(([a], [b]) => a.localeCompare(b)).map(([subject, hit]) => (
                                  <li
                                    key={subject}
                                    className="pl-14 pr-4 py-1.5 text-xs break-all"
                                  >
                                    <div>
                                      {isUniprotAccession(subject)
                                        ? <UniprotProteinLabel rawId={subject} />
                                        : <span className="font-mono text-emerald-700 dark:text-emerald-400">{subject}</span>
                                      }
                                    </div>
                                    {(hit.evalue !== null || hit.pident !== null) && (
                                      <div className="mt-0.5 flex gap-3 text-gray-400 dark:text-gray-500">
                                        {hit.evalue !== null && (
                                          <span>E: <span className="font-mono text-gray-500 dark:text-gray-400">{hit.evalue.toExponential(1)}</span></span>
                                        )}
                                        {hit.pident !== null && (
                                          <span>ID: <span className="font-mono text-gray-500 dark:text-gray-400">{hit.pident.toFixed(1)}%</span></span>
                                        )}
                                      </div>
                                    )}
                                  </li>
                                ))}
                              </ul>
                            )}
                          </li>
                        )
                      })}
                    </ul>
                  )}
                </li>
              )
            })}
          </ul>
        )}
      </div>

      {hasSelection && !loading && hierarchy.size > 0 && (
        <div className="px-4 py-2 border-t border-gray-100 dark:border-gray-700 bg-gray-50 dark:bg-gray-800 text-xs text-gray-500 dark:text-gray-400">
          {hierarchy.size} peptide{hierarchy.size !== 1 ? 's' : ''}
        </div>
      )}
    </div>
  )
}
