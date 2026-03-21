import { describe, it, expect } from 'vitest'
import { parseComboCsv, comboRowsToTaxonNodes, comboRowsToGoTermNodes } from '../comboParser'

const HEADER = 'tax_id,tax_name,tax_rank,parent_tax_id,go_id,go_name,go_namespace,parent_go_ids,quantity,fraction_of_taxon,fraction_of_go,ratio_total_taxon,ratio_total_go,n_peptides'
const ENRICHED_HEADER = `${HEADER},pvalue_go_for_taxon,pvalue_taxon_for_go,qvalue_go_for_taxon,qvalue_taxon_for_go,zscore_go_for_taxon,zscore_taxon_for_go`

describe('parseComboCsv', () => {
  it('returns empty array for empty input', () => {
    expect(parseComboCsv('')).toEqual([])
  })

  it('returns empty array for header-only input', () => {
    expect(parseComboCsv(HEADER)).toEqual([])
  })

  it('parses a single row', () => {
    const csv = `${HEADER}\n30,ClassA,class,20,GO:0000004,C,biological_process,GO:0000002;GO:0000003,10.0,0.5,0.8,0.12,0.04,3`
    const rows = parseComboCsv(csv)
    expect(rows).toHaveLength(1)
    expect(rows[0]).toEqual({
      taxId: '30',
      taxName: 'ClassA',
      taxRank: 'class',
      parentTaxId: '20',
      goId: 'GO:0000004',
      goName: 'C',
      goNamespace: 'biological_process',
      parentGoIds: ['GO:0000002', 'GO:0000003'],
      quantity: 10.0,
      fractionOfTaxon: 0.5,
      fractionOfGo: 0.8,
      ratioTotalTaxon: 0.12,
      ratioTotalGo: 0.04,
      nPeptides: 3,
    })
  })

  it('parses multiple rows', () => {
    const csv = [
      HEADER,
      '30,ClassA,class,20,GO:0000004,C,biological_process,GO:0000002;GO:0000003,10.0,0.5,0.8,0.12,0.04,3',
      '20,PhylumA,phylum,10,GO:0000002,A,biological_process,GO:0000001,20.0,1.0,0.6,0.20,0.08,5',
    ].join('\n')
    const rows = parseComboCsv(csv)
    expect(rows).toHaveLength(2)
  })

  it('handles empty parent_go_ids', () => {
    const csv = `${HEADER}\n30,ClassA,class,20,GO:0000001,root,biological_process,,10.0,0.5,0.8,0.12,0.04,3`
    const rows = parseComboCsv(csv)
    expect(rows[0].parentGoIds).toEqual([])
  })

  it('handles empty parent_tax_id (root)', () => {
    const csv = `${HEADER}\n1,root,root,,GO:0000001,root_go,biological_process,,10.0,1.0,1.0,1.0,1.0,1`
    const rows = parseComboCsv(csv)
    expect(rows[0].parentTaxId).toBe('')
  })

  it('skips rows with fewer than 14 fields', () => {
    const csv = [
      HEADER,
      '30,ClassA,class,20,GO:0000004,C,biological_process,GO:0000002;GO:0000003,10.0,0.5,0.8,0.12,0.04,3',
      '30,ClassA,class,20,short',
    ].join('\n')
    const rows = parseComboCsv(csv)
    expect(rows).toHaveLength(1)
  })

  it('parses enrichment columns when present', () => {
    const csv = `${ENRICHED_HEADER}\n30,ClassA,class,20,GO:0000004,C,biological_process,GO:0000002;GO:0000003,10.0,0.5,0.8,0.12,0.04,3,0.01,0.02,0.03,0.04,1.5,-2.5`
    const rows = parseComboCsv(csv)
    expect(rows[0].pvalueGoForTaxon).toBe(0.01)
    expect(rows[0].pvalueTaxonForGo).toBe(0.02)
    expect(rows[0].qvalueGoForTaxon).toBe(0.03)
    expect(rows[0].qvalueTaxonForGo).toBe(0.04)
    expect(rows[0].zscoreGoForTaxon).toBe(1.5)
    expect(rows[0].zscoreTaxonForGo).toBe(-2.5)
  })

  it('parses signed infinite z-scores when present', () => {
    const csv = `${ENRICHED_HEADER}\n30,ClassA,class,20,GO:0000004,C,biological_process,GO:0000002;GO:0000003,10.0,0.5,0.8,0.12,0.04,3,0.01,0.02,0.03,0.04,+inf,-inf`
    const rows = parseComboCsv(csv)
    expect(rows[0].zscoreGoForTaxon).toBe(Number.POSITIVE_INFINITY)
    expect(rows[0].zscoreTaxonForGo).toBe(Number.NEGATIVE_INFINITY)
  })

  it('keeps enrichment columns undefined for old CSVs', () => {
    const csv = `${HEADER}\n30,ClassA,class,20,GO:0000004,C,biological_process,GO:0000002;GO:0000003,10.0,0.5,0.8,0.12,0.04,3`
    const rows = parseComboCsv(csv)
    expect(rows[0].qvalueGoForTaxon).toBeUndefined()
    expect(rows[0].zscoreTaxonForGo).toBeUndefined()
  })
})

describe('comboRowsToTaxonNodes', () => {
  const rows = parseComboCsv([
    HEADER,
    '30,ClassA,class,20,GO:0000004,C,biological_process,GO:0000002;GO:0000003,10.0,0.5,0.8,0.12,0.04,3',
    '20,PhylumA,phylum,10,GO:0000004,C,biological_process,GO:0000002;GO:0000003,20.0,1.0,0.6,0.20,0.04,5',
    '30,ClassA,class,20,GO:0000002,A,biological_process,GO:0000001,5.0,0.25,0.3,0.12,0.08,2',
  ].join('\n'))

  it('filters to the specified GO term', () => {
    const nodes = comboRowsToTaxonNodes(rows, 'GO:0000004')
    expect(nodes).toHaveLength(2)
    expect(nodes.map(n => n.taxId).sort()).toEqual(['20', '30'])
  })

  it('maps fields correctly', () => {
    const nodes = comboRowsToTaxonNodes(parseComboCsv([
      ENRICHED_HEADER,
      '30,ClassA,class,20,GO:0000004,C,biological_process,GO:0000002;GO:0000003,10.0,0.5,0.8,0.12,0.04,3,0.01,0.02,0.03,0.04,1.5,-2.5',
    ].join('\n')), 'GO:0000004')
    const node30 = nodes.find(n => n.taxId === '30')!
    expect(node30.name).toBe('ClassA')
    expect(node30.rank).toBe('class')
    expect(node30.parentTaxId).toBe('20')
    expect(node30.quantity).toBe(10.0)
    expect(node30.ratioTotal).toBe(0.12)
    expect(node30.ratioAnnotated).toBe(0)
    expect(node30.nPeptides).toBe(3)
    expect(node30.fractionOfTaxon).toBe(0.5)
    expect(node30.fractionOfGo).toBe(0.8)
    expect(node30.qvalueTaxonForGo).toBe(0.04)
    expect(node30.zscoreTaxonForGo).toBe(-2.5)
  })

  it('returns empty for non-existent GO term', () => {
    const nodes = comboRowsToTaxonNodes(rows, 'GO:9999999')
    expect(nodes).toHaveLength(0)
  })
})

describe('comboRowsToGoTermNodes', () => {
  const rows = parseComboCsv([
    HEADER,
    '30,ClassA,class,20,GO:0000004,C,biological_process,GO:0000002;GO:0000003,10.0,0.5,0.8,0.12,0.04,3',
    '30,ClassA,class,20,GO:0000002,A,biological_process,GO:0000001,5.0,0.25,0.3,0.12,0.08,2',
    '20,PhylumA,phylum,10,GO:0000004,C,biological_process,GO:0000002;GO:0000003,20.0,1.0,0.6,0.20,0.04,5',
  ].join('\n'))

  it('filters to the specified tax ID', () => {
    const nodes = comboRowsToGoTermNodes(rows, '30')
    expect(nodes).toHaveLength(2)
    expect(nodes.map(n => n.id).sort()).toEqual(['GO:0000002', 'GO:0000004'])
  })

  it('maps fields correctly', () => {
    const nodes = comboRowsToGoTermNodes(parseComboCsv([
      ENRICHED_HEADER,
      '30,ClassA,class,20,GO:0000004,C,biological_process,GO:0000002;GO:0000003,10.0,0.5,0.8,0.12,0.04,3,0.01,0.02,0.03,0.04,1.5,-2.5',
    ].join('\n')), '30')
    const node = nodes.find(n => n.id === 'GO:0000004')!
    expect(node.name).toBe('C')
    expect(node.namespace).toBe('biological_process')
    expect(node.parentIds).toEqual(['GO:0000002', 'GO:0000003'])
    expect(node.quantity).toBe(10.0)
    expect(node.ratioTotal).toBe(0.04)
    expect(node.ratioAnnotated).toBe(0)
    expect(node.nPeptides).toBe(3)
    expect(node.fractionOfTaxon).toBe(0.5)
    expect(node.fractionOfGo).toBe(0.8)
    expect(node.qvalueGoForTaxon).toBe(0.03)
    expect(node.zscoreGoForTaxon).toBe(1.5)
  })

  it('returns empty for non-existent tax ID', () => {
    const nodes = comboRowsToGoTermNodes(rows, '999')
    expect(nodes).toHaveLength(0)
  })
})
