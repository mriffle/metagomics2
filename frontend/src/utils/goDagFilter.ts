import type { GoTermNode } from './goParser'

/**
 * Iteratively prune insignificant *leaf* GO terms from a DAG by q-value.
 *
 * A node is a "leaf" when no other present node lists it as a parent. Leaves whose
 * `qvalueGoForTaxon` is worse than `maxQvalue` (i.e. `q > maxQvalue`) are removed.
 * Removing a leaf can expose its parent as a new leaf, which is then re-evaluated,
 * so this repeats to a fixpoint: every remaining leaf passes the threshold.
 *
 * Internal nodes are never pruned by q-value. A node with a passing descendant
 * always retains a child and so is preserved, even if its own q-value is worse than
 * the threshold. This keeps the DAG spine from significant leaves up to the root(s)
 * intact. When every node along a branch fails, the branch is removed entirely
 * (including a root that becomes an isolated failing leaf) — the logical end of
 * "prune leaves until all remaining leaves pass".
 *
 * A leaf whose q-value is missing (`null`/`undefined`) is treated as not prunable
 * and kept. In practice this cannot happen while the filter is usable — see the
 * `showQvalueMetric` gate in GoDagPage, and note that q-value eligibility is an
 * all-or-nothing property of the selected taxon — but the guard keeps the function
 * total.
 *
 * The relative order of surviving nodes is preserved.
 */
export function pruneInsignificantLeaves(
  nodes: GoTermNode[],
  maxQvalue: number,
): GoTermNode[] {
  const present = new Set(nodes.map((n) => n.id))
  const byId = new Map(nodes.map((n) => [n.id, n] as const))

  // Distinct, present parents per node, and childCount[id] = number of present
  // nodes that list `id` as a parent (i.e. `id`'s in-graph child count).
  const parentsOf = new Map<string, string[]>()
  const childCount = new Map<string, number>()
  for (const id of present) childCount.set(id, 0)
  for (const node of nodes) {
    const parents = [...new Set(node.parentIds)].filter((pid) => present.has(pid))
    parentsOf.set(node.id, parents)
    for (const pid of parents) {
      childCount.set(pid, (childCount.get(pid) ?? 0) + 1)
    }
  }

  // Worklist seeded with the initial leaves (no present children).
  const queue: string[] = []
  for (const node of nodes) {
    if ((childCount.get(node.id) ?? 0) === 0) queue.push(node.id)
  }

  const removed = new Set<string>()
  while (queue.length > 0) {
    const id = queue.pop()!
    if (removed.has(id) || (childCount.get(id) ?? 0) !== 0) continue // not (still) a leaf
    const q = byId.get(id)!.qvalueGoForTaxon
    // Keep leaves that pass, or that have no q-value to judge.
    if (q == null || q <= maxQvalue) continue
    // Prune this failing leaf; a parent becomes a new leaf once all its children go.
    removed.add(id)
    for (const pid of parentsOf.get(id) ?? []) {
      if (removed.has(pid)) continue
      const next = (childCount.get(pid) ?? 0) - 1
      childCount.set(pid, next)
      if (next === 0) queue.push(pid)
    }
  }

  if (removed.size === 0) return nodes
  return nodes.filter((n) => !removed.has(n.id))
}
