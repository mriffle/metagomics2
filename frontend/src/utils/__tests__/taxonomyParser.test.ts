import { describe, it, expect } from 'vitest'
import { parseTaxonomyCsv, filterCanonicalRanks, filterByMaxRank, validateCanonicalHierarchy, ensureStrictRankLayers, CANONICAL_RANKS, CANONICAL_RANKS_ORDERED } from '../taxonomyParser'
import type { TaxonNode } from '../taxonomyParser'

const HEADER = 'tax_id,name,rank,parent_tax_id,quantity,ratio_total,ratio_annotated,n_peptides'

describe('parseTaxonomyCsv', () => {
  it('returns empty array for empty input', () => {
    expect(parseTaxonomyCsv('')).toEqual([])
  })

  it('returns empty array for header-only input', () => {
    expect(parseTaxonomyCsv(HEADER)).toEqual([])
  })

  it('parses a single row', () => {
    const csv = `${HEADER}\n1,root,root,,1000.0,1.0,1.0,100`
    const nodes = parseTaxonomyCsv(csv)
    expect(nodes).toHaveLength(1)
    expect(nodes[0]).toEqual({
      taxId: '1',
      name: 'root',
      rank: 'root',
      parentTaxId: '',
      quantity: 1000.0,
      ratioTotal: 1.0,
      ratioAnnotated: 1.0,
      nPeptides: 100,
    })
  })

  it('parses multiple rows', () => {
    const csv = [
      HEADER,
      '1,root,root,,1000.0,1.0,1.0,100',
      '2,Bacteria,domain,1,800.0,0.8,0.9,80',
      '9606,Homo sapiens,species,9605,50.0,0.05,0.1,5',
    ].join('\n')
    const nodes = parseTaxonomyCsv(csv)
    expect(nodes).toHaveLength(3)
    expect(nodes[2].taxId).toBe('9606')
    expect(nodes[2].name).toBe('Homo sapiens')
    expect(nodes[2].parentTaxId).toBe('9605')
  })

  it('handles empty parent_tax_id for root', () => {
    const csv = `${HEADER}\n1,root,root,,1000.0,1.0,1.0,100`
    const nodes = parseTaxonomyCsv(csv)
    expect(nodes[0].parentTaxId).toBe('')
  })

  it('handles high-precision decimal values', () => {
    const csv = `${HEADER}\n1,root,root,,1234.5678900000,0.0500000000,0.1000000000,42`
    const nodes = parseTaxonomyCsv(csv)
    expect(nodes[0].quantity).toBeCloseTo(1234.56789, 5)
    expect(nodes[0].ratioTotal).toBeCloseTo(0.05, 10)
  })

  it('skips rows with fewer than 8 fields', () => {
    const csv = [
      HEADER,
      '1,root,root,,1000.0,1.0,1.0,100',
      '2,short,row',
    ].join('\n')
    const nodes = parseTaxonomyCsv(csv)
    expect(nodes).toHaveLength(1)
  })

  it('handles quoted names with commas', () => {
    const csv = `${HEADER}\n9606,"Homo sapiens, neanderthalensis",subspecies,9605,50.0,0.05,0.1,5`
    const nodes = parseTaxonomyCsv(csv)
    expect(nodes[0].name).toBe('Homo sapiens, neanderthalensis')
  })
})

describe('CANONICAL_RANKS', () => {
  it('contains all expected ranks', () => {
    const expected = ['root', 'domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    for (const rank of expected) {
      expect(CANONICAL_RANKS.has(rank)).toBe(true)
    }
  })

  it('does not contain non-canonical ranks', () => {
    expect(CANONICAL_RANKS.has('no rank')).toBe(false)
    expect(CANONICAL_RANKS.has('clade')).toBe(false)
    expect(CANONICAL_RANKS.has('tribe')).toBe(false)
    expect(CANONICAL_RANKS.has('subphylum')).toBe(false)
    expect(CANONICAL_RANKS.has('superfamily')).toBe(false)
  })
})

describe('filterCanonicalRanks', () => {
  function makeNode(taxId: string, name: string, rank: string, parentTaxId: string): TaxonNode {
    return { taxId, name, rank, parentTaxId, quantity: 10, ratioTotal: 0.1, ratioAnnotated: 0.2, nPeptides: 1 }
  }

  it('keeps only canonical rank nodes', () => {
    const nodes: TaxonNode[] = [
      makeNode('1', 'root', 'root', ''),
      makeNode('2', 'Bacteria', 'domain', '1'),
      makeNode('3', 'SomeClade', 'clade', '2'),
      makeNode('4', 'Proteobacteria', 'phylum', '3'),
    ]
    const result = filterCanonicalRanks(nodes)
    expect(result).toHaveLength(3)
    expect(result.map(n => n.taxId)).toEqual(['1', '2', '4'])
  })

  it('re-links parent through non-canonical intermediate', () => {
    const nodes: TaxonNode[] = [
      makeNode('1', 'root', 'root', ''),
      makeNode('2', 'Bacteria', 'domain', '1'),
      makeNode('3', 'SomeClade', 'clade', '2'),
      makeNode('4', 'Proteobacteria', 'phylum', '3'),
    ]
    const result = filterCanonicalRanks(nodes)
    // Phylum (4) should link to domain (2), skipping clade (3)
    const phylum = result.find(n => n.taxId === '4')!
    expect(phylum.parentTaxId).toBe('2')
  })

  it('re-links through multiple non-canonical intermediates', () => {
    const nodes: TaxonNode[] = [
      makeNode('1', 'root', 'root', ''),
      makeNode('2', 'Bacteria', 'domain', '1'),
      makeNode('3', 'clade1', 'clade', '2'),
      makeNode('4', 'clade2', 'no rank', '3'),
      makeNode('5', 'clade3', 'clade', '4'),
      makeNode('6', 'Proteobacteria', 'phylum', '5'),
    ]
    const result = filterCanonicalRanks(nodes)
    const phylum = result.find(n => n.taxId === '6')!
    expect(phylum.parentTaxId).toBe('2')
  })

  it('preserves parent when already canonical', () => {
    const nodes: TaxonNode[] = [
      makeNode('1', 'root', 'root', ''),
      makeNode('2', 'Bacteria', 'domain', '1'),
      makeNode('3', 'Proteobacteria', 'phylum', '2'),
    ]
    const result = filterCanonicalRanks(nodes)
    const phylum = result.find(n => n.taxId === '3')!
    expect(phylum.parentTaxId).toBe('2')
  })

  it('handles root with empty parent', () => {
    const nodes: TaxonNode[] = [
      makeNode('1', 'root', 'root', ''),
      makeNode('2', 'Bacteria', 'domain', '1'),
    ]
    const result = filterCanonicalRanks(nodes)
    const root = result.find(n => n.taxId === '1')!
    expect(root.parentTaxId).toBe('')
  })

  it('handles parent not in data (dangling reference)', () => {
    const nodes: TaxonNode[] = [
      makeNode('100', 'SomeSpecies', 'species', '99'),
    ]
    const result = filterCanonicalRanks(nodes)
    expect(result).toHaveLength(1)
    // Parent 99 not in data, should become empty
    expect(result[0].parentTaxId).toBe('')
  })

  it('handles non-canonical parent chain ending at missing node', () => {
    const nodes: TaxonNode[] = [
      makeNode('1', 'root', 'root', ''),
      makeNode('3', 'SomeClade', 'clade', '2'),  // parent 2 not in data
      makeNode('4', 'Proteobacteria', 'phylum', '3'),
    ]
    const result = filterCanonicalRanks(nodes)
    const phylum = result.find(n => n.taxId === '4')!
    // Walking up: 3 is clade (non-canonical), parent is 2 which is missing → empty
    expect(phylum.parentTaxId).toBe('')
  })

  it('does not modify original nodes', () => {
    const nodes: TaxonNode[] = [
      makeNode('1', 'root', 'root', ''),
      makeNode('2', 'clade', 'clade', '1'),
      makeNode('3', 'Bacteria', 'domain', '2'),
    ]
    const originalParent = nodes[2].parentTaxId
    filterCanonicalRanks(nodes)
    expect(nodes[2].parentTaxId).toBe(originalParent)
  })

  it('returns empty array when no canonical nodes exist', () => {
    const nodes: TaxonNode[] = [
      makeNode('1', 'clade1', 'clade', ''),
      makeNode('2', 'clade2', 'no rank', '1'),
    ]
    const result = filterCanonicalRanks(nodes)
    expect(result).toEqual([])
  })

  it('normalizes NCBI root (no rank) to canonical root', () => {
    const nodes: TaxonNode[] = [
      makeNode('1', 'root', 'no rank', ''),
      makeNode('131567', 'cellular organisms', 'no rank', '1'),
      makeNode('2', 'Bacteria', 'domain', '131567'),
    ]
    const result = filterCanonicalRanks(nodes)
    expect(result).toHaveLength(2)

    const root = result.find(n => n.taxId === '1')!
    expect(root.rank).toBe('root')

    const bacteria = result.find(n => n.taxId === '2')!
    expect(bacteria.parentTaxId).toBe('1')
  })

  it('handles a realistic taxonomy hierarchy', () => {
    const nodes: TaxonNode[] = [
      makeNode('1', 'root', 'no rank', ''),
      makeNode('131567', 'cellular organisms', 'no rank', '1'),
      makeNode('2', 'Bacteria', 'domain', '131567'),
      makeNode('1224', 'Pseudomonadota', 'phylum', '2'),
      makeNode('1236', 'Gammaproteobacteria', 'class', '1224'),
      makeNode('72274', 'unclassified Gammaproteobacteria', 'no rank', '1236'),
      makeNode('91347', 'Enterobacterales', 'order', '1236'),
      makeNode('543', 'Enterobacteriaceae', 'family', '91347'),
      makeNode('561', 'Escherichia', 'genus', '543'),
      makeNode('562', 'Escherichia coli', 'species', '561'),
    ]
    const result = filterCanonicalRanks(nodes)

    // Should have: root, Bacteria, Pseudomonadota, Gammaproteobacteria,
    // Enterobacterales, Enterobacteriaceae, Escherichia, E. coli
    expect(result).toHaveLength(8)

    // Bacteria should link to root (skipping 'cellular organisms')
    const bacteria = result.find(n => n.taxId === '2')!
    expect(bacteria.parentTaxId).toBe('1')

    // E. coli should link to Escherichia
    const ecoli = result.find(n => n.taxId === '562')!
    expect(ecoli.parentTaxId).toBe('561')

    // Enterobacterales should link to Gammaproteobacteria
    const enterobacterales = result.find(n => n.taxId === '91347')!
    expect(enterobacterales.parentTaxId).toBe('1236')
  })
})

describe('CANONICAL_RANKS_ORDERED', () => {
  it('is ordered from root to species', () => {
    expect(CANONICAL_RANKS_ORDERED[0]).toBe('root')
    expect(CANONICAL_RANKS_ORDERED[CANONICAL_RANKS_ORDERED.length - 1]).toBe('species')
  })

  it('has the same elements as CANONICAL_RANKS set', () => {
    expect(CANONICAL_RANKS_ORDERED.length).toBe(CANONICAL_RANKS.size)
    for (const r of CANONICAL_RANKS_ORDERED) {
      expect(CANONICAL_RANKS.has(r)).toBe(true)
    }
  })
})

describe('filterByMaxRank', () => {
  function makeNode(taxId: string, name: string, rank: string, parentTaxId: string): TaxonNode {
    return { taxId, name, rank, parentTaxId, quantity: 10, ratioTotal: 0.1, ratioAnnotated: 0.2, nPeptides: 1 }
  }

  const fullHierarchy: TaxonNode[] = [
    makeNode('1', 'root', 'root', ''),
    makeNode('2', 'Bacteria', 'domain', '1'),
    makeNode('3', 'Fungi', 'kingdom', '1'),
    makeNode('4', 'Proteobacteria', 'phylum', '2'),
    makeNode('5', 'Gammaproteobacteria', 'class', '4'),
    makeNode('6', 'Enterobacterales', 'order', '5'),
    makeNode('7', 'Enterobacteriaceae', 'family', '6'),
    makeNode('8', 'Escherichia', 'genus', '7'),
    makeNode('9', 'Escherichia coli', 'species', '8'),
  ]

  it('filters to domain (root + domain)', () => {
    const result = filterByMaxRank(fullHierarchy, 'domain')
    expect(result.map(n => n.rank)).toEqual(['root', 'domain'])
  })

  it('filters to phylum (root + domain + kingdom + phylum)', () => {
    const result = filterByMaxRank(fullHierarchy, 'phylum')
    const ranks = result.map(n => n.rank)
    expect(ranks).toEqual(['root', 'domain', 'kingdom', 'phylum'])
  })

  it('filters to class', () => {
    const result = filterByMaxRank(fullHierarchy, 'class')
    const ranks = result.map(n => n.rank)
    expect(ranks).toEqual(['root', 'domain', 'kingdom', 'phylum', 'class'])
  })

  it('filters to species (all nodes)', () => {
    const result = filterByMaxRank(fullHierarchy, 'species')
    expect(result.length).toBe(fullHierarchy.length)
  })

  it('filters to root only', () => {
    const result = filterByMaxRank(fullHierarchy, 'root')
    expect(result).toHaveLength(1)
    expect(result[0].rank).toBe('root')
  })

  it('preserves node data unchanged', () => {
    const result = filterByMaxRank(fullHierarchy, 'phylum')
    const phylum = result.find(n => n.taxId === '4')!
    expect(phylum.name).toBe('Proteobacteria')
    expect(phylum.parentTaxId).toBe('2')
    expect(phylum.quantity).toBe(10)
  })

  it('returns empty array when no nodes match', () => {
    const result = filterByMaxRank([], 'species')
    expect(result).toEqual([])
  })
})

describe('validateCanonicalHierarchy', () => {
  function makeNode(taxId: string, name: string, rank: string, parentTaxId: string): TaxonNode {
    return { taxId, name, rank, parentTaxId, quantity: 10, ratioTotal: 0.1, ratioAnnotated: 0.2, nPeptides: 1 }
  }

  it('returns no errors for a valid consecutive hierarchy', () => {
    const nodes: TaxonNode[] = [
      makeNode('1', 'root', 'root', ''),
      makeNode('2', 'Bacteria', 'domain', '1'),
      makeNode('3', 'Proteobacteria', 'kingdom', '2'),
      makeNode('4', 'Gammaproteobacteria', 'phylum', '3'),
    ]
    expect(validateCanonicalHierarchy(nodes)).toEqual([])
  })

  it('allows rank gaps (phylum parent is domain, skipping kingdom)', () => {
    const nodes: TaxonNode[] = [
      makeNode('1', 'root', 'root', ''),
      makeNode('2', 'Bacteria', 'domain', '1'),
      makeNode('4', 'Proteobacteria', 'phylum', '2'),
    ]
    expect(validateCanonicalHierarchy(nodes)).toEqual([])
  })

  it('allows rank gaps (phylum parent is root, skipping domain+kingdom)', () => {
    const nodes: TaxonNode[] = [
      makeNode('1', 'root', 'root', ''),
      makeNode('4', 'Proteobacteria', 'phylum', '1'),
    ]
    expect(validateCanonicalHierarchy(nodes)).toEqual([])
  })

  it('returns error when parent is missing from data', () => {
    const nodes: TaxonNode[] = [
      makeNode('1', 'root', 'root', ''),
      makeNode('4', 'Proteobacteria', 'phylum', '99'),
    ]
    const errors = validateCanonicalHierarchy(nodes)
    expect(errors.length).toBeGreaterThan(0)
    expect(errors[0]).toContain('99')
    expect(errors[0]).toContain('not in the data')
  })

  it('returns error when non-root node has no parent', () => {
    const nodes: TaxonNode[] = [
      makeNode('1', 'root', 'root', ''),
      makeNode('2', 'Bacteria', 'domain', ''),
    ]
    const errors = validateCanonicalHierarchy(nodes)
    expect(errors.length).toBeGreaterThan(0)
    expect(errors[0]).toContain('has no parent')
  })

  it('returns error when parent is not a higher rank', () => {
    const nodes: TaxonNode[] = [
      makeNode('1', 'root', 'root', ''),
      makeNode('2', 'Bacteria', 'domain', '1'),
      makeNode('3', 'SomePhylum', 'phylum', '2'),
      makeNode('4', 'BadDomain', 'domain', '3'),
    ]
    const errors = validateCanonicalHierarchy(nodes)
    expect(errors.length).toBe(1)
    expect(errors[0]).toContain('not a higher rank')
  })

  it('returns error when parent is same rank', () => {
    const nodes: TaxonNode[] = [
      makeNode('1', 'root', 'root', ''),
      makeNode('2', 'Bacteria', 'domain', '1'),
      makeNode('3', 'Archaea', 'domain', '2'),
    ]
    const errors = validateCanonicalHierarchy(nodes)
    expect(errors.length).toBe(1)
    expect(errors[0]).toContain('not a higher rank')
  })

  it('does not report errors for root node', () => {
    const nodes: TaxonNode[] = [
      makeNode('1', 'root', 'root', ''),
    ]
    expect(validateCanonicalHierarchy(nodes)).toEqual([])
  })

  it('returns no errors for empty input', () => {
    expect(validateCanonicalHierarchy([])).toEqual([])
  })

  it('validates a complete realistic hierarchy with no errors', () => {
    const nodes: TaxonNode[] = [
      makeNode('1', 'root', 'root', ''),
      makeNode('2', 'Bacteria', 'domain', '1'),
      makeNode('3', 'SomeKingdom', 'kingdom', '2'),
      makeNode('4', 'Proteobacteria', 'phylum', '3'),
      makeNode('5', 'Gammaproteobacteria', 'class', '4'),
      makeNode('6', 'Enterobacterales', 'order', '5'),
      makeNode('7', 'Enterobacteriaceae', 'family', '6'),
      makeNode('8', 'Escherichia', 'genus', '7'),
      makeNode('9', 'Escherichia coli', 'species', '8'),
    ]
    expect(validateCanonicalHierarchy(nodes)).toEqual([])
  })
})

describe('ensureStrictRankLayers', () => {
  function makeNode(taxId: string, name: string, rank: string, parentTaxId: string, quantity = 10): TaxonNode {
    return { taxId, name, rank, parentTaxId, quantity, ratioTotal: 0.1, ratioAnnotated: 0.2, nPeptides: 1 }
  }

  it('does not modify nodes when all ranks are consecutive', () => {
    const nodes: TaxonNode[] = [
      makeNode('1', 'root', 'root', ''),
      makeNode('2', 'Bacteria', 'domain', '1'),
      makeNode('3', 'Fungi', 'kingdom', '2'),
    ]
    const result = ensureStrictRankLayers(nodes)
    expect(result).toHaveLength(3)
    expect(result.filter(n => n.taxId.startsWith('__placeholder'))).toHaveLength(0)
  })

  it('inserts placeholder when phylum parent is domain (skipping kingdom)', () => {
    const nodes: TaxonNode[] = [
      makeNode('1', 'root', 'root', ''),
      makeNode('2', 'Bacteria', 'domain', '1'),
      makeNode('4', 'Proteobacteria', 'phylum', '2'),
    ]
    const result = ensureStrictRankLayers(nodes)
    expect(result).toHaveLength(4)
    const ranks = result.map(n => n.rank)
    expect(ranks).toEqual(['root', 'domain', 'kingdom', 'phylum'])

    const phylum = result.find(n => n.taxId === '4')!
    expect(phylum.parentTaxId).not.toBe('2')

    const placeholder = result.find(n => n.rank === 'kingdom' && n.taxId.startsWith('__placeholder'))!
    expect(placeholder.parentTaxId).toBe('2')
    expect(placeholder.name).toBe('(no kingdom for Proteobacteria)')
    expect(phylum.parentTaxId).toBe(placeholder.taxId)
  })

  it('inserts multiple placeholders for large gaps', () => {
    const nodes: TaxonNode[] = [
      makeNode('1', 'root', 'root', ''),
      makeNode('2', 'Bacteria', 'domain', '1'),
      makeNode('5', 'Gammaproteobacteria', 'class', '2'),
    ]
    const result = ensureStrictRankLayers(nodes)
    const ranks = result.map(n => n.rank)
    expect(ranks).toEqual(['root', 'domain', 'kingdom', 'phylum', 'class'])
  })

  it('handles multiple children with the same gap', () => {
    const nodes: TaxonNode[] = [
      makeNode('1', 'root', 'root', ''),
      makeNode('2', 'Eukaryota', 'domain', '1'),
      makeNode('4', 'Haptophyta', 'phylum', '2', 50),
      makeNode('5', 'Rhodophyta', 'phylum', '2', 30),
    ]
    const result = ensureStrictRankLayers(nodes)
    const phyla = result.filter(n => n.rank === 'phylum' && !n.taxId.startsWith('__placeholder'))
    expect(phyla).toHaveLength(2)

    // Each phylum should be at depth 3 (root=0, domain=1, kingdom=2, phylum=3)
    const nodeMap = new Map(result.map(n => [n.taxId, n]))
    for (const p of phyla) {
      let depth = 0
      let current = p
      while (current.parentTaxId) {
        const par = nodeMap.get(current.parentTaxId)
        if (!par) break
        current = par
        depth++
      }
      expect(depth).toBe(3)
    }
  })

  it('returns empty array for empty input', () => {
    expect(ensureStrictRankLayers([])).toEqual([])
  })

  it('placeholder nodes have descriptive name', () => {
    const nodes: TaxonNode[] = [
      makeNode('1', 'root', 'root', ''),
      makeNode('2', 'Bacteria', 'domain', '1'),
      makeNode('4', 'Proteobacteria', 'phylum', '2'),
    ]
    const result = ensureStrictRankLayers(nodes)
    const placeholders = result.filter(n => n.taxId.startsWith('__placeholder'))
    expect(placeholders).toHaveLength(1)
    expect(placeholders[0].name).toBe('(no kingdom for Proteobacteria)')
  })
})
