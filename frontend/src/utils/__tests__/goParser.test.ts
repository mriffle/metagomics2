import { describe, it, expect } from 'vitest'
import { parseGoTermsCsv } from '../goParser'

const HEADER = 'go_id,name,namespace,parent_go_ids,quantity,ratio_total,ratio_annotated,n_peptides'

describe('parseGoTermsCsv', () => {
  it('returns empty array for empty input', () => {
    expect(parseGoTermsCsv('')).toEqual([])
  })

  it('returns empty array for header-only input', () => {
    expect(parseGoTermsCsv(HEADER)).toEqual([])
  })

  it('parses a single row', () => {
    const csv = `${HEADER}\nGO:0000001,root_BP,biological_process,,100.0,0.5,1.0,10`
    const nodes = parseGoTermsCsv(csv)
    expect(nodes).toHaveLength(1)
    expect(nodes[0]).toEqual({
      id: 'GO:0000001',
      name: 'root_BP',
      namespace: 'biological_process',
      parentIds: [],
      quantity: 100.0,
      ratioTotal: 0.5,
      ratioAnnotated: 1.0,
      nPeptides: 10,
    })
  })

  it('parses multiple rows', () => {
    const csv = [
      HEADER,
      'GO:0000001,root_BP,biological_process,,100.0,0.5,1.0,10',
      'GO:0000002,A,biological_process,GO:0000001,50.0,0.25,0.5,5',
    ].join('\n')
    const nodes = parseGoTermsCsv(csv)
    expect(nodes).toHaveLength(2)
    expect(nodes[1].id).toBe('GO:0000002')
    expect(nodes[1].parentIds).toEqual(['GO:0000001'])
  })

  it('parses semicolon-delimited parent IDs', () => {
    const csv = `${HEADER}\nGO:0000004,C,biological_process,GO:0000002;GO:0000003,10.0,0.1,0.2,5`
    const nodes = parseGoTermsCsv(csv)
    expect(nodes[0].parentIds).toEqual(['GO:0000002', 'GO:0000003'])
  })

  it('handles empty parent_go_ids', () => {
    const csv = `${HEADER}\nGO:0000001,root,biological_process,,100.0,0.5,1.0,10`
    const nodes = parseGoTermsCsv(csv)
    expect(nodes[0].parentIds).toEqual([])
  })

  it('trims whitespace from fields', () => {
    const csv = `${HEADER}\n GO:0000001 , root_BP , biological_process , , 100.0 , 0.5 , 1.0 , 10`
    const nodes = parseGoTermsCsv(csv)
    expect(nodes[0].id).toBe('GO:0000001')
    expect(nodes[0].name).toBe('root_BP')
    expect(nodes[0].namespace).toBe('biological_process')
  })

  it('handles quoted names with commas', () => {
    const csv = `${HEADER}\nGO:0000001,"regulation of something, important",biological_process,,100.0,0.5,1.0,10`
    const nodes = parseGoTermsCsv(csv)
    expect(nodes[0].name).toBe('regulation of something, important')
  })

  it('defaults to 0 for non-numeric quantity', () => {
    const csv = `${HEADER}\nGO:0000001,root,biological_process,,bad,0.5,1.0,10`
    const nodes = parseGoTermsCsv(csv)
    expect(nodes[0].quantity).toBe(0)
  })

  it('defaults to 0 for non-numeric nPeptides', () => {
    const csv = `${HEADER}\nGO:0000001,root,biological_process,,100.0,0.5,1.0,bad`
    const nodes = parseGoTermsCsv(csv)
    expect(nodes[0].nPeptides).toBe(0)
  })

  it('skips rows with fewer than 8 fields', () => {
    const csv = [
      HEADER,
      'GO:0000001,root,biological_process,,100.0,0.5,1.0,10',
      'GO:0000002,short,row',
    ].join('\n')
    const nodes = parseGoTermsCsv(csv)
    expect(nodes).toHaveLength(1)
  })

  it('handles high-precision decimal values', () => {
    const csv = `${HEADER}\nGO:0000001,root,biological_process,,1234.5678900000,0.0500000000,0.1000000000,42`
    const nodes = parseGoTermsCsv(csv)
    expect(nodes[0].quantity).toBeCloseTo(1234.56789, 5)
    expect(nodes[0].ratioTotal).toBeCloseTo(0.05, 10)
    expect(nodes[0].ratioAnnotated).toBeCloseTo(0.1, 10)
  })

  it('handles Windows-style line endings', () => {
    const csv = `${HEADER}\r\nGO:0000001,root,biological_process,,100.0,0.5,1.0,10\r\n`
    const nodes = parseGoTermsCsv(csv)
    expect(nodes).toHaveLength(1)
    expect(nodes[0].id).toBe('GO:0000001')
  })
})
