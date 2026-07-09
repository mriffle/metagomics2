import { describe, it, expect } from 'vitest'
import { pruneInsignificantLeaves } from '../goDagFilter'
import type { GoTermNode } from '../goParser'

function node(id: string, parentIds: string[], qvalue?: number): GoTermNode {
  return {
    id,
    name: id,
    namespace: 'biological_process',
    parentIds,
    quantity: 0,
    ratioTotal: 0,
    ratioAnnotated: 0,
    nPeptides: 0,
    qvalueGoForTaxon: qvalue,
  }
}

function ids(nodes: GoTermNode[]): string[] {
  return nodes.map((n) => n.id)
}

describe('pruneInsignificantLeaves', () => {
  it('returns an empty array for empty input', () => {
    expect(pruneInsignificantLeaves([], 0.05)).toEqual([])
  })

  it('keeps everything when all leaves pass', () => {
    // root -> A -> B, both A and B pass (root has no q-value but is never a leaf)
    const nodes = [
      node('root', []),
      node('A', ['root'], 0.01),
      node('B', ['A'], 0.02),
    ]
    expect(ids(pruneInsignificantLeaves(nodes, 0.05))).toEqual(['root', 'A', 'B'])
  })

  it('removes a failing leaf and re-checks the exposed parent', () => {
    // root -> A(pass) -> B(fail): B removed, A becomes a passing leaf and is kept.
    const nodes = [
      node('root', []),
      node('A', ['root'], 0.01),
      node('B', ['A'], 0.5),
    ]
    expect(ids(pruneInsignificantLeaves(nodes, 0.05))).toEqual(['root', 'A'])
  })

  it('cascades all the way to the root when every node on a branch fails', () => {
    const nodes = [
      node('root', [], 0.9),
      node('A', ['root'], 0.8),
      node('B', ['A'], 0.7),
    ]
    expect(pruneInsignificantLeaves(nodes, 0.05)).toEqual([])
  })

  it('preserves ancestors with bad q-values when a descendant passes', () => {
    // root(bad) -> A(bad) -> B(good): B is a passing leaf, so A and root are kept
    // despite their own q-values being worse than the threshold.
    const nodes = [
      node('root', [], 0.9),
      node('A', ['root'], 0.8),
      node('B', ['A'], 0.001),
    ]
    expect(ids(pruneInsignificantLeaves(nodes, 0.05))).toEqual(['root', 'A', 'B'])
  })

  it('exposes both parents in a diamond when the shared leaf is pruned', () => {
    // root -> {A, B} -> C. C is the only leaf; pruning it exposes A and B.
    // A fails and is removed; B passes and is kept, keeping root alive too.
    const nodes = [
      node('root', []),
      node('A', ['root'], 0.9),
      node('B', ['root'], 0.01),
      node('C', ['A', 'B'], 0.9),
    ]
    expect(ids(pruneInsignificantLeaves(nodes, 0.05))).toEqual(['root', 'B'])
  })

  it('keeps a parent that still has a passing sibling leaf', () => {
    // P -> C1(pass), P -> C2(fail): C2 removed, P retains child C1 and is kept.
    const nodes = [
      node('P', []),
      node('C1', ['P'], 0.01),
      node('C2', ['P'], 0.9),
    ]
    expect(ids(pruneInsignificantLeaves(nodes, 0.05))).toEqual(['P', 'C1'])
  })

  it('treats the threshold as inclusive (q == max passes)', () => {
    const nodes = [node('root', []), node('A', ['root'], 0.05)]
    expect(ids(pruneInsignificantLeaves(nodes, 0.05))).toEqual(['root', 'A'])
  })

  it('prunes a leaf just over the threshold', () => {
    const nodes = [node('root', []), node('A', ['root'], 0.0500001)]
    expect(ids(pruneInsignificantLeaves(nodes, 0.05))).toEqual(['root'])
  })

  it('keeps a leaf whose q-value is missing (unjudgeable)', () => {
    const nodes = [
      node('root', []),
      node('A', ['root'], undefined),
      node('B', ['root'], 0.9),
    ]
    // A has no q-value → kept; B fails → pruned; root retains child A → kept.
    expect(ids(pruneInsignificantLeaves(nodes, 0.05))).toEqual(['root', 'A'])
  })

  it('prunes an isolated failing node (both root and leaf)', () => {
    expect(pruneInsignificantLeaves([node('X', [], 0.9)], 0.05)).toEqual([])
  })

  it('keeps an isolated passing node', () => {
    expect(ids(pruneInsignificantLeaves([node('X', [], 0.01)], 0.05))).toEqual(['X'])
  })

  it('ignores parent references to absent nodes (dangling parents)', () => {
    // A's parent "ghost" is not present; A is still a leaf and root-like.
    const nodes = [node('A', ['ghost'], 0.9)]
    expect(pruneInsignificantLeaves(nodes, 0.05)).toEqual([])
  })

  it('preserves input order among survivors', () => {
    const nodes = [
      node('root', []),
      node('A', ['root'], 0.9), // leaf, fails
      node('B', ['root'], 0.01), // leaf, passes
      node('C', ['root'], 0.02), // leaf, passes
    ]
    expect(ids(pruneInsignificantLeaves(nodes, 0.05))).toEqual(['root', 'B', 'C'])
  })

  it('returns the original array reference when nothing is pruned', () => {
    const nodes = [node('root', []), node('A', ['root'], 0.01)]
    expect(pruneInsignificantLeaves(nodes, 0.05)).toBe(nodes)
  })
})
