import { parseCSVLine } from './csvParser'

export interface TaxonNode {
  taxId: string
  name: string
  rank: string
  parentTaxId: string
  quantity: number
  ratioTotal: number
  ratioAnnotated: number
  nPeptides: number
}

export const CANONICAL_RANKS_ORDERED = [
  'root',
  'domain',
  'kingdom',
  'phylum',
  'class',
  'order',
  'family',
  'genus',
  'species',
] as const

export type CanonicalRank = typeof CANONICAL_RANKS_ORDERED[number]

export const CANONICAL_RANKS = new Set<string>(CANONICAL_RANKS_ORDERED)

export function parseTaxonomyCsv(text: string): TaxonNode[] {
  const lines = text.trim().split('\n')
  if (lines.length < 2) return []

  const nodes: TaxonNode[] = []
  for (let i = 1; i < lines.length; i++) {
    const fields = parseCSVLine(lines[i])
    if (fields.length < 8) continue

    const [taxId, name, rank, parentTaxId, quantity, ratioTotal, ratioAnnotated, nPeptides] = fields

    nodes.push({
      taxId: taxId.trim(),
      name: name.trim(),
      rank: rank.trim(),
      parentTaxId: parentTaxId.trim(),
      quantity: parseFloat(quantity) || 0,
      ratioTotal: parseFloat(ratioTotal) || 0,
      ratioAnnotated: parseFloat(ratioAnnotated) || 0,
      nPeptides: parseInt(nPeptides, 10) || 0,
    })
  }
  return nodes
}

/**
 * Filter canonical nodes to only include ranks up to (and including) maxRank.
 */
export function filterByMaxRank(canonicalNodes: TaxonNode[], maxRank: CanonicalRank): TaxonNode[] {
  const maxIdx = CANONICAL_RANKS_ORDERED.indexOf(maxRank)
  const allowedRanks = new Set<string>(CANONICAL_RANKS_ORDERED.slice(0, maxIdx + 1))

  return canonicalNodes.filter(n => allowedRanks.has(n.rank))
}

/**
 * Build a rank index map for O(1) lookups.
 */
const RANK_INDEX = new Map<string, number>()
for (let i = 0; i < CANONICAL_RANKS_ORDERED.length; i++) {
  RANK_INDEX.set(CANONICAL_RANKS_ORDERED[i], i)
}

/**
 * Validate that the canonical hierarchy is well-formed.
 *
 * Checks:
 * 1. Every non-root node's parentTaxId references another node in the set.
 * 2. Every non-root node's parent is at a strictly higher (lower index)
 *    canonical rank. Rank gaps are allowed (NCBI doesn't always have every
 *    canonical rank in every lineage).
 *
 * Returns an array of error messages (empty if valid).
 */
export function validateCanonicalHierarchy(nodes: TaxonNode[]): string[] {
  const errors: string[] = []
  const nodeMap = new Map<string, TaxonNode>()
  for (const n of nodes) nodeMap.set(n.taxId, n)

  for (const node of nodes) {
    if (node.rank === 'root') continue

    // Check 1: parent must exist in the set
    if (!node.parentTaxId) {
      errors.push(`${node.name} (${node.taxId}, rank=${node.rank}) has no parent`)
      continue
    }

    const parent = nodeMap.get(node.parentTaxId)
    if (!parent) {
      errors.push(`${node.name} (${node.taxId}, rank=${node.rank}) references parent ${node.parentTaxId} which is not in the data`)
      continue
    }

    // Check 2: parent must be at a strictly higher rank
    const nodeIdx = RANK_INDEX.get(node.rank)
    const parentIdx = RANK_INDEX.get(parent.rank)
    if (nodeIdx !== undefined && parentIdx !== undefined) {
      if (parentIdx >= nodeIdx) {
        errors.push(
          `${node.name} (${node.taxId}, rank=${node.rank}) has parent ${parent.name} (${parent.taxId}, rank=${parent.rank}) which is not a higher rank`
        )
      }
    }
  }

  return errors
}

/**
 * Ensure strict rank layering for sunburst/treemap charts.
 *
 * Plotly determines ring depth by the parent-child chain length, not by any
 * rank metadata. If a phylum's canonical parent is a domain (because there's
 * no kingdom in that lineage), it would appear in ring 3 instead of ring 4.
 *
 * This function inserts invisible placeholder nodes for every missing
 * intermediate canonical rank between a node and its parent, so the
 * parent-child chain depth always matches the rank depth.
 *
 * Placeholder nodes have:
 *  - taxId: `__placeholder_<rank>_<childTaxId>`
 *  - name: empty string (invisible in the chart)
 *  - quantity/ratios: same as the child (so branchvalues='total' works)
 */
export function ensureStrictRankLayers(nodes: TaxonNode[]): TaxonNode[] {
  const nodeMap = new Map<string, TaxonNode>()
  for (const n of nodes) nodeMap.set(n.taxId, n)

  const result: TaxonNode[] = []
  const added = new Set<string>()

  for (const node of nodes) {
    const nodeRankIdx = RANK_INDEX.get(node.rank)
    if (nodeRankIdx === undefined) { result.push(node); continue }

    const parent = node.parentTaxId ? nodeMap.get(node.parentTaxId) : undefined
    const parentRankIdx = parent ? RANK_INDEX.get(parent.rank) : undefined

    // Determine the effective parent rank index (-1 if no parent = above root)
    const effectiveParentIdx = parentRankIdx !== undefined ? parentRankIdx : -1
    const effectiveParentId = parentRankIdx !== undefined ? node.parentTaxId : ''

    // If parent is exactly one rank above (or root with no gap), no placeholders needed
    if (nodeRankIdx - effectiveParentIdx <= 1) {
      result.push(node)
      continue
    }

    // Insert placeholders for each missing intermediate rank
    let currentParentId = effectiveParentId
    for (let ri = effectiveParentIdx + 1; ri < nodeRankIdx; ri++) {
      const placeholderRank = CANONICAL_RANKS_ORDERED[ri]
      const placeholderId = `__placeholder_${placeholderRank}_${node.taxId}`

      if (!added.has(placeholderId)) {
        added.add(placeholderId)
        result.push({
          taxId: placeholderId,
          name: `(no ${placeholderRank} for ${node.name})`,
          rank: placeholderRank,
          parentTaxId: currentParentId,
          quantity: node.quantity,
          ratioTotal: node.ratioTotal,
          ratioAnnotated: node.ratioAnnotated,
          nPeptides: node.nPeptides,
        })
      }
      currentParentId = placeholderId
    }

    // Re-parent the actual node to the last placeholder
    result.push({
      ...node,
      parentTaxId: currentParentId,
    })
  }

  return result
}

/**
 * Filter to canonical ranks and re-link parents.
 * For each canonical node, walk up the parent chain until finding
 * the nearest ancestor that is also at a canonical rank.
 */
export function filterCanonicalRanks(allNodes: TaxonNode[]): TaxonNode[] {
  const nodeMap = new Map<string, TaxonNode>()
  for (const n of allNodes) nodeMap.set(n.taxId, n)

  // Find the true root node (no parent or self-referencing parent) and
  // normalize its rank to 'root' so it's treated as canonical regardless
  // of its actual rank string (NCBI uses "no rank" for the root).
  // Only include it if there are other canonical nodes that could use it.
  const hasCanonicalNodes = allNodes.some(n => CANONICAL_RANKS.has(n.rank))
  const rootNode = hasCanonicalNodes
    ? allNodes.find(n => !n.parentTaxId || n.parentTaxId === n.taxId)
    : undefined
  const canonicalNodes = allNodes
    .filter(n => CANONICAL_RANKS.has(n.rank) || (rootNode && n.taxId === rootNode.taxId))
    .map(n => (rootNode && n.taxId === rootNode.taxId && n.rank !== 'root')
      ? { ...n, rank: 'root' }
      : n
    )

  // For each canonical node, find its nearest canonical ancestor
  return canonicalNodes.map(node => {
    let parentId = node.parentTaxId
    const visited = new Set<string>()

    while (parentId && !visited.has(parentId)) {
      visited.add(parentId)
      const parent = nodeMap.get(parentId)
      if (!parent) {
        parentId = ''
        break
      }
      if (CANONICAL_RANKS.has(parent.rank) || (rootNode && parent.taxId === rootNode.taxId)) {
        break
      }
      parentId = parent.parentTaxId
    }

    return {
      ...node,
      parentTaxId: parentId,
    }
  })
}
