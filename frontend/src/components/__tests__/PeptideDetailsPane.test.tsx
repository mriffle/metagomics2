import { describe, it, expect, vi, beforeEach } from 'vitest'
import { render, screen, fireEvent, waitFor } from '@testing-library/react'
import PeptideDetailsPane from '../PeptideDetailsPane'

const { mockQuery } = vi.hoisted(() => ({ mockQuery: vi.fn() }))

vi.mock('../../utils/duckdb', () => ({
  getDuckDB: vi.fn().mockResolvedValue({
    db: {},
    conn: { query: mockQuery },
  }),
  registerMappingFile: vi.fn().mockResolvedValue(undefined),
}))

function makeQueryResult(rows: Record<string, unknown>[]) {
  return {
    toArray: () =>
      rows.map((row) => ({
        toJSON: () => row,
      })),
  }
}

describe('PeptideDetailsPane', () => {
  beforeEach(() => {
    vi.clearAllMocks()
    mockQuery.mockResolvedValue(makeQueryResult([]))
  })

  it('renders placeholder when no selectedTaxId or selectedGoId', () => {
    render(<PeptideDetailsPane jobId="job1" listId="list_000" />)
    expect(screen.getByText(/click a node/i)).toBeTruthy()
  })

  it('renders placeholder when both selectedTaxIds and selectedGoId are null', () => {
    render(
      <PeptideDetailsPane
        jobId="job1"
        listId="list_000"
        selectedTaxIds={null}
        selectedGoId={null}
      />
    )
    expect(screen.getByText(/click a node/i)).toBeTruthy()
  })

  it('shows loading state while DuckDB query is in flight', async () => {
    let resolveQuery!: (value: unknown) => void
    mockQuery.mockReturnValue(new Promise((resolve) => { resolveQuery = resolve }))

    render(
      <PeptideDetailsPane jobId="job1" listId="list_000" selectedTaxIds={['9606']} />
    )

    await waitFor(() => {
      expect(screen.getByText(/loading/i)).toBeTruthy()
    })

    resolveQuery(makeQueryResult([]))
  })

  it('renders correct peptide hierarchy when query returns rows', async () => {
    mockQuery.mockResolvedValue(
      makeQueryResult([
        {
          peptide: 'PEPTIDESEQ',
          peptide_lca_tax_ids: [9606, 1],
          peptide_go_terms: ['GO:0000001'],
          background_protein: 'prot_A',
          annotated_protein: 'subj_X',
          evalue: 1e-20,
          pident: 95.5,
        },
      ])
    )

    render(
      <PeptideDetailsPane jobId="job1" listId="list_000" selectedTaxIds={['9606']} />
    )

    await waitFor(() => {
      expect(screen.getByText('PEPTIDESEQ')).toBeTruthy()
    })
  })

  it('shows "no peptides found" message when query returns empty results', async () => {
    mockQuery.mockResolvedValue(makeQueryResult([]))

    render(
      <PeptideDetailsPane jobId="job1" listId="list_000" selectedTaxIds={['9606']} />
    )

    await waitFor(() => {
      expect(screen.getByText(/no peptides found/i)).toBeTruthy()
    })
  })

  it('expands peptide row to show background proteins on click', async () => {
    mockQuery.mockResolvedValue(
      makeQueryResult([
        {
          peptide: 'PEPTIDESEQ',
          peptide_lca_tax_ids: [9606, 1],
          peptide_go_terms: [],
          background_protein: 'prot_A',
          annotated_protein: 'subj_X',
          evalue: 1e-10,
          pident: 90.0,
        },
      ])
    )

    render(
      <PeptideDetailsPane jobId="job1" listId="list_000" selectedTaxIds={['9606']} />
    )

    await waitFor(() => { expect(screen.getByText('PEPTIDESEQ')).toBeTruthy() })

    fireEvent.click(screen.getByText('PEPTIDESEQ').closest('button')!)

    await waitFor(() => {
      expect(screen.getByText('prot_A')).toBeTruthy()
    })
  })

  it('collapses peptide row on second click', async () => {
    mockQuery.mockResolvedValue(
      makeQueryResult([
        {
          peptide: 'PEPTIDESEQ',
          peptide_lca_tax_ids: [9606, 1],
          peptide_go_terms: [],
          background_protein: 'prot_A',
          annotated_protein: 'subj_X',
          evalue: 1e-10,
          pident: 90.0,
        },
      ])
    )

    render(
      <PeptideDetailsPane jobId="job1" listId="list_000" selectedTaxIds={['9606']} />
    )

    await waitFor(() => { expect(screen.getByText('PEPTIDESEQ')).toBeTruthy() })

    const peptideButton = screen.getByText('PEPTIDESEQ').closest('button')!
    fireEvent.click(peptideButton)
    await waitFor(() => { expect(screen.getByText('prot_A')).toBeTruthy() })

    fireEvent.click(peptideButton)
    await waitFor(() => {
      expect(screen.queryByText('prot_A')).toBeNull()
    })
  })

  it('expands background protein to show annotated proteins on click', async () => {
    mockQuery.mockResolvedValue(
      makeQueryResult([
        {
          peptide: 'PEPTIDESEQ',
          peptide_lca_tax_ids: [9606, 1],
          peptide_go_terms: [],
          background_protein: 'prot_A',
          annotated_protein: 'subj_X',
          evalue: 1e-10,
          pident: 90.0,
        },
      ])
    )

    render(
      <PeptideDetailsPane jobId="job1" listId="list_000" selectedTaxIds={['9606']} />
    )

    await waitFor(() => { expect(screen.getByText('PEPTIDESEQ')).toBeTruthy() })
    fireEvent.click(screen.getByText('PEPTIDESEQ').closest('button')!)

    await waitFor(() => { expect(screen.getByText('prot_A')).toBeTruthy() })
    fireEvent.click(screen.getByText('prot_A').closest('button')!)

    await waitFor(() => {
      expect(screen.getByText('subj_X')).toBeTruthy()
    })
  })

  it('shows selection info in header when taxIds are provided', async () => {
    mockQuery.mockResolvedValue(makeQueryResult([]))

    render(
      <PeptideDetailsPane jobId="job1" listId="list_000" selectedTaxIds={['9606']} />
    )

    await waitFor(() => {
      expect(screen.getByText(/Tax ID: 9606/)).toBeTruthy()
    })
  })

  it('shows GO ID in header when goId is provided', async () => {
    mockQuery.mockResolvedValue(makeQueryResult([]))

    render(
      <PeptideDetailsPane jobId="job1" listId="list_000" selectedGoId="GO:0000001" />
    )

    await waitFor(() => {
      expect(screen.getByText('GO:0000001')).toBeTruthy()
    })
  })
})
