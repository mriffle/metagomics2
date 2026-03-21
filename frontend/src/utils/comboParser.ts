import { parseCSVLine } from './csvParser'
import type { TaxonNode } from './taxonomyParser'
import type { GoTermNode } from './goParser'

export interface ComboRow {
  taxId: string
  taxName: string
  taxRank: string
  parentTaxId: string
  goId: string
  goName: string
  goNamespace: string
  parentGoIds: string[]
  quantity: number
  fractionOfTaxon: number
  fractionOfGo: number
  ratioTotalTaxon: number
  ratioTotalGo: number
  nPeptides: number
  pvalueGoForTaxon?: number
  pvalueTaxonForGo?: number
  qvalueGoForTaxon?: number
  qvalueTaxonForGo?: number
  zscoreGoForTaxon?: number
  zscoreTaxonForGo?: number
}

/**
 * Parse go_taxonomy_combo.csv into ComboRow[].
 */
export function parseComboCsv(text: string): ComboRow[] {
  const lines = text.trim().split('\n')
  if (lines.length < 2) return []

  const headerFields = parseCSVLine(lines[0]).map(field => field.trim())
  const headerIndex = new Map<string, number>()
  headerFields.forEach((field, index) => headerIndex.set(field, index))

  const requiredHeaders = [
    'tax_id',
    'tax_name',
    'tax_rank',
    'parent_tax_id',
    'go_id',
    'go_name',
    'go_namespace',
    'parent_go_ids',
    'quantity',
    'fraction_of_taxon',
    'fraction_of_go',
    'ratio_total_taxon',
    'ratio_total_go',
    'n_peptides',
  ]
  if (requiredHeaders.some(header => !headerIndex.has(header))) return []

  const getField = (fields: string[], header: string): string => {
    const index = headerIndex.get(header)
    return index == null ? '' : (fields[index] ?? '')
  }

  const parseOptionalNumber = (fields: string[], header: string): number | undefined => {
    const raw = getField(fields, header).trim()
    if (!raw) return undefined
    if (/^\+?inf(?:inity)?$/i.test(raw)) return Number.POSITIVE_INFINITY
    if (/^-inf(?:inity)?$/i.test(raw)) return Number.NEGATIVE_INFINITY
    const parsed = parseFloat(raw)
    return Number.isNaN(parsed) ? undefined : parsed
  }

  const rows: ComboRow[] = []
  for (let i = 1; i < lines.length; i++) {
    const fields = parseCSVLine(lines[i])
    if (fields.length < requiredHeaders.length) continue

    const taxId = getField(fields, 'tax_id')
    const taxName = getField(fields, 'tax_name')
    const taxRank = getField(fields, 'tax_rank')
    const parentTaxId = getField(fields, 'parent_tax_id')
    const goId = getField(fields, 'go_id')
    const goName = getField(fields, 'go_name')
    const goNamespace = getField(fields, 'go_namespace')
    const parentGoIds = getField(fields, 'parent_go_ids')
    const quantity = getField(fields, 'quantity')
    const fractionOfTaxon = getField(fields, 'fraction_of_taxon')
    const fractionOfGo = getField(fields, 'fraction_of_go')
    const ratioTotalTaxon = getField(fields, 'ratio_total_taxon')
    const ratioTotalGo = getField(fields, 'ratio_total_go')
    const nPeptides = getField(fields, 'n_peptides')

    rows.push({
      taxId: taxId.trim(),
      taxName: taxName.trim(),
      taxRank: taxRank.trim(),
      parentTaxId: parentTaxId.trim(),
      goId: goId.trim(),
      goName: goName.trim(),
      goNamespace: goNamespace.trim(),
      parentGoIds: parentGoIds.trim()
        ? parentGoIds.trim().split(';').map(s => s.trim()).filter(Boolean)
        : [],
      quantity: parseFloat(quantity) || 0,
      fractionOfTaxon: parseFloat(fractionOfTaxon) || 0,
      fractionOfGo: parseFloat(fractionOfGo) || 0,
      ratioTotalTaxon: parseFloat(ratioTotalTaxon) || 0,
      ratioTotalGo: parseFloat(ratioTotalGo) || 0,
      nPeptides: parseInt(nPeptides, 10) || 0,
      pvalueGoForTaxon: parseOptionalNumber(fields, 'pvalue_go_for_taxon'),
      pvalueTaxonForGo: parseOptionalNumber(fields, 'pvalue_taxon_for_go'),
      qvalueGoForTaxon: parseOptionalNumber(fields, 'qvalue_go_for_taxon'),
      qvalueTaxonForGo: parseOptionalNumber(fields, 'qvalue_taxon_for_go'),
      zscoreGoForTaxon: parseOptionalNumber(fields, 'zscore_go_for_taxon'),
      zscoreTaxonForGo: parseOptionalNumber(fields, 'zscore_taxon_for_go'),
    })
  }
  return rows
}

/**
 * Filter combo rows by a GO term and reshape into TaxonNode[].
 *
 * Uses `ratioTotalTaxon` as `ratioTotal` so the min-ratio filter still
 * works against the taxon's fraction of total quantity.
 * `fractionOfTaxon` is stored separately for tooltip display.
 */
export function comboRowsToTaxonNodes(rows: ComboRow[], goId: string): TaxonNode[] {
  const byTaxId = new Map<string, ComboRow>()
  for (const row of rows) {
    if (row.goId !== goId) continue
    byTaxId.set(row.taxId, row)
  }

  return Array.from(byTaxId.values()).map(row => ({
    taxId: row.taxId,
    name: row.taxName,
    rank: row.taxRank,
    parentTaxId: row.parentTaxId,
    quantity: row.quantity,
    ratioTotal: row.ratioTotalTaxon,
    ratioAnnotated: 0,
    nPeptides: row.nPeptides,
    fractionOfTaxon: row.fractionOfTaxon,
    fractionOfGo: row.fractionOfGo,
    qvalueTaxonForGo: row.qvalueTaxonForGo,
    zscoreTaxonForGo: row.zscoreTaxonForGo,
  }))
}

/**
 * Filter combo rows by a taxonomy node and reshape into GoTermNode[].
 *
 * Uses `ratioTotalGo` as `ratioTotal` so the min-ratio filter still
 * works against the GO term's fraction of total quantity.
 * `fractionOfGo` is stored separately for tooltip display.
 */
export function comboRowsToGoTermNodes(rows: ComboRow[], taxId: string): GoTermNode[] {
  const byGoId = new Map<string, ComboRow>()
  for (const row of rows) {
    if (row.taxId !== taxId) continue
    byGoId.set(row.goId, row)
  }

  return Array.from(byGoId.values()).map(row => ({
    id: row.goId,
    name: row.goName,
    namespace: row.goNamespace,
    parentIds: row.parentGoIds,
    quantity: row.quantity,
    ratioTotal: row.ratioTotalGo,
    ratioAnnotated: 0,
    nPeptides: row.nPeptides,
    fractionOfTaxon: row.fractionOfTaxon,
    fractionOfGo: row.fractionOfGo,
    qvalueGoForTaxon: row.qvalueGoForTaxon,
    zscoreGoForTaxon: row.zscoreGoForTaxon,
  }))
}
