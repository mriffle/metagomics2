import { describe, it, expect } from 'vitest'
import { parseCSVLine } from '../csvParser'

describe('parseCSVLine', () => {
  it('parses simple comma-separated fields', () => {
    expect(parseCSVLine('a,b,c')).toEqual(['a', 'b', 'c'])
  })

  it('parses single field', () => {
    expect(parseCSVLine('hello')).toEqual(['hello'])
  })

  it('parses empty string as single empty field', () => {
    expect(parseCSVLine('')).toEqual([''])
  })

  it('handles empty fields between commas', () => {
    expect(parseCSVLine('a,,c')).toEqual(['a', '', 'c'])
  })

  it('handles trailing comma', () => {
    expect(parseCSVLine('a,b,')).toEqual(['a', 'b', ''])
  })

  it('handles leading comma', () => {
    expect(parseCSVLine(',a,b')).toEqual(['', 'a', 'b'])
  })

  it('parses quoted fields', () => {
    expect(parseCSVLine('"hello","world"')).toEqual(['hello', 'world'])
  })

  it('handles commas inside quoted fields', () => {
    expect(parseCSVLine('"hello, world",foo')).toEqual(['hello, world', 'foo'])
  })

  it('handles escaped quotes (doubled)', () => {
    expect(parseCSVLine('"say ""hello""",bar')).toEqual(['say "hello"', 'bar'])
  })

  it('handles mixed quoted and unquoted fields', () => {
    expect(parseCSVLine('plain,"quoted",another')).toEqual(['plain', 'quoted', 'another'])
  })

  it('handles newlines inside quoted fields', () => {
    expect(parseCSVLine('"line1\nline2",b')).toEqual(['line1\nline2', 'b'])
  })

  it('parses a realistic taxonomy CSV line', () => {
    const line = '9606,Homo sapiens,species,9605,1234.5678900000,0.0500000000,0.1000000000,42'
    const fields = parseCSVLine(line)
    expect(fields).toEqual([
      '9606', 'Homo sapiens', 'species', '9605',
      '1234.5678900000', '0.0500000000', '0.1000000000', '42',
    ])
  })

  it('parses a realistic GO CSV line with semicolon-delimited parents', () => {
    const line = 'GO:0000004,"cell morphogenesis",biological_process,GO:0000002;GO:0000003,10.0000000000,0.1000000000,0.2000000000,5'
    const fields = parseCSVLine(line)
    expect(fields[0]).toBe('GO:0000004')
    expect(fields[1]).toBe('cell morphogenesis')
    expect(fields[3]).toBe('GO:0000002;GO:0000003')
  })

  it('handles quoted field with name containing comma', () => {
    const line = '9606,"Homo sapiens, neanderthalensis",subspecies,9605,100.0,0.01,0.02,3'
    const fields = parseCSVLine(line)
    expect(fields[1]).toBe('Homo sapiens, neanderthalensis')
    expect(fields.length).toBe(8)
  })
})
