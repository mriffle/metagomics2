import { parseCSVLine } from './csvParser'

export interface GoTermNode {
  id: string
  name: string
  namespace: string
  parentIds: string[]
  quantity: number
  ratioTotal: number
  ratioAnnotated: number
  nPeptides: number
}

export function parseGoTermsCsv(text: string): GoTermNode[] {
  const lines = text.trim().split('\n')
  if (lines.length < 2) return []

  const nodes: GoTermNode[] = []
  for (let i = 1; i < lines.length; i++) {
    const fields = parseCSVLine(lines[i])
    if (fields.length < 8) continue

    const [goId, name, namespace, parentGoIds, quantity, ratioTotal, ratioAnnotated, nPeptides] = fields

    nodes.push({
      id: goId.trim(),
      name: name.trim(),
      namespace: namespace.trim(),
      parentIds: parentGoIds.trim()
        ? parentGoIds.trim().split(';').map(s => s.trim()).filter(Boolean)
        : [],
      quantity: parseFloat(quantity) || 0,
      ratioTotal: parseFloat(ratioTotal) || 0,
      ratioAnnotated: parseFloat(ratioAnnotated) || 0,
      nPeptides: parseInt(nPeptides, 10) || 0,
    })
  }
  return nodes
}
