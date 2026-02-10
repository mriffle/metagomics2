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
}

/**
 * Parse go_taxonomy_combo.csv into ComboRow[].
 */
export function parseComboCsv(text: string): ComboRow[] {
  const lines = text.trim().split('\n')
  if (lines.length < 2) return []

  const rows: ComboRow[] = []
  for (let i = 1; i < lines.length; i++) {
    const fields = parseCSVLine(lines[i])
    if (fields.length < 14) continue

    const [
      taxId, taxName, taxRank, parentTaxId,
      goId, goName, goNamespace, parentGoIds,
      quantity, fractionOfTaxon, fractionOfGo,
      ratioTotalTaxon, ratioTotalGo, nPeptides,
    ] = fields

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
  }))
}
