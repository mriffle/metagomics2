import { useState } from 'react'
import { useNavigate } from 'react-router-dom'
import { Upload, FileText, Settings, Loader2, AlertCircle } from 'lucide-react'

export default function NewJobPage() {
  const navigate = useNavigate()
  const [fastaFile, setFastaFile] = useState<File | null>(null)
  const [peptideFiles, setPeptideFiles] = useState<File[]>([])
  const [submitting, setSubmitting] = useState(false)
  const [error, setError] = useState<string | null>(null)

  // Parameters
  const [searchTool, setSearchTool] = useState('diamond')
  const [maxEvalue, setMaxEvalue] = useState('')
  const [minPident, setMinPident] = useState('')
  const [topK, setTopK] = useState('')

  async function handleSubmit(e: React.FormEvent) {
    e.preventDefault()

    if (!fastaFile) {
      setError('Please select a FASTA file')
      return
    }
    if (peptideFiles.length === 0) {
      setError('Please select at least one peptide file')
      return
    }

    setSubmitting(true)
    setError(null)

    try {
      const formData = new FormData()
      formData.append('fasta', fastaFile)
      peptideFiles.forEach((file) => {
        formData.append('peptides', file)
      })

      const params: Record<string, unknown> = {
        search_tool: searchTool,
      }
      if (maxEvalue) params.max_evalue = parseFloat(maxEvalue)
      if (minPident) params.min_pident = parseFloat(minPident)
      if (topK) params.top_k = parseInt(topK)

      formData.append('params', JSON.stringify(params))

      const response = await fetch('/api/jobs', {
        method: 'POST',
        body: formData,
      })

      if (!response.ok) {
        const data = await response.json()
        throw new Error(data.detail || 'Failed to create job')
      }

      const data = await response.json()
      navigate(`/job/${data.job_id}`)
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Unknown error')
    } finally {
      setSubmitting(false)
    }
  }

  return (
    <div className="max-w-2xl mx-auto">
      <div className="mb-8">
        <h1 className="text-3xl font-bold text-gray-900">New Annotation Job</h1>
        <p className="mt-2 text-gray-600">
          Upload your data and configure parameters for annotation
        </p>
      </div>

      {error && (
        <div className="mb-6 p-4 bg-red-50 border border-red-200 rounded-lg flex items-center gap-2 text-red-700">
          <AlertCircle className="w-5 h-5 flex-shrink-0" />
          {error}
        </div>
      )}

      <form onSubmit={handleSubmit} className="space-y-8">
        {/* File Uploads */}
        <div className="bg-white rounded-lg border border-gray-200 p-6">
          <h2 className="text-lg font-semibold text-gray-900 mb-4 flex items-center gap-2">
            <Upload className="w-5 h-5" />
            Input Files
          </h2>

          <div className="space-y-4">
            {/* FASTA Upload */}
            <div>
              <label className="block text-sm font-medium text-gray-700 mb-2">
                Background Proteome FASTA *
              </label>
              <div className="border-2 border-dashed border-gray-300 rounded-lg p-4 text-center hover:border-indigo-400 transition-colors">
                <input
                  type="file"
                  accept=".fasta,.fa,.faa"
                  onChange={(e) => setFastaFile(e.target.files?.[0] || null)}
                  className="hidden"
                  id="fasta-upload"
                />
                <label htmlFor="fasta-upload" className="cursor-pointer">
                  {fastaFile ? (
                    <div className="flex items-center justify-center gap-2 text-indigo-600">
                      <FileText className="w-5 h-5" />
                      {fastaFile.name}
                    </div>
                  ) : (
                    <div className="text-gray-500">
                      <Upload className="w-8 h-8 mx-auto mb-2" />
                      <span>Click to upload FASTA file</span>
                    </div>
                  )}
                </label>
              </div>
            </div>

            {/* Peptide Files Upload */}
            <div>
              <label className="block text-sm font-medium text-gray-700 mb-2">
                Peptide Lists (CSV/TSV) *
              </label>
              <div className="border-2 border-dashed border-gray-300 rounded-lg p-4 text-center hover:border-indigo-400 transition-colors">
                <input
                  type="file"
                  accept=".csv,.tsv,.txt"
                  multiple
                  onChange={(e) => setPeptideFiles(Array.from(e.target.files || []))}
                  className="hidden"
                  id="peptide-upload"
                />
                <label htmlFor="peptide-upload" className="cursor-pointer">
                  {peptideFiles.length > 0 ? (
                    <div className="text-indigo-600">
                      <FileText className="w-5 h-5 mx-auto mb-2" />
                      {peptideFiles.length} file(s) selected
                      <ul className="text-sm mt-2">
                        {peptideFiles.map((f, i) => (
                          <li key={i}>{f.name}</li>
                        ))}
                      </ul>
                    </div>
                  ) : (
                    <div className="text-gray-500">
                      <Upload className="w-8 h-8 mx-auto mb-2" />
                      <span>Click to upload peptide files</span>
                    </div>
                  )}
                </label>
              </div>
            </div>
          </div>
        </div>

        {/* Parameters */}
        <div className="bg-white rounded-lg border border-gray-200 p-6">
          <h2 className="text-lg font-semibold text-gray-900 mb-4 flex items-center gap-2">
            <Settings className="w-5 h-5" />
            Parameters
          </h2>

          <div className="grid grid-cols-2 gap-4">
            <div>
              <label className="block text-sm font-medium text-gray-700 mb-1">
                Search Tool
              </label>
              <select
                value={searchTool}
                onChange={(e) => setSearchTool(e.target.value)}
                className="w-full px-3 py-2 border border-gray-300 rounded-lg focus:ring-2 focus:ring-indigo-500 focus:border-indigo-500"
              >
                <option value="diamond">DIAMOND</option>
                <option value="blast">BLAST</option>
              </select>
            </div>

            <div>
              <label className="block text-sm font-medium text-gray-700 mb-1">
                Max E-value
              </label>
              <input
                type="text"
                value={maxEvalue}
                onChange={(e) => setMaxEvalue(e.target.value)}
                placeholder="e.g., 1e-5"
                className="w-full px-3 py-2 border border-gray-300 rounded-lg focus:ring-2 focus:ring-indigo-500 focus:border-indigo-500"
              />
            </div>

            <div>
              <label className="block text-sm font-medium text-gray-700 mb-1">
                Min % Identity
              </label>
              <input
                type="text"
                value={minPident}
                onChange={(e) => setMinPident(e.target.value)}
                placeholder="e.g., 80"
                className="w-full px-3 py-2 border border-gray-300 rounded-lg focus:ring-2 focus:ring-indigo-500 focus:border-indigo-500"
              />
            </div>

            <div>
              <label className="block text-sm font-medium text-gray-700 mb-1">
                Top K Hits
              </label>
              <input
                type="text"
                value={topK}
                onChange={(e) => setTopK(e.target.value)}
                placeholder="e.g., 10"
                className="w-full px-3 py-2 border border-gray-300 rounded-lg focus:ring-2 focus:ring-indigo-500 focus:border-indigo-500"
              />
            </div>
          </div>
        </div>

        {/* Submit */}
        <div className="flex justify-end">
          <button
            type="submit"
            disabled={submitting}
            className="inline-flex items-center gap-2 px-6 py-3 bg-indigo-600 text-white rounded-lg hover:bg-indigo-700 disabled:opacity-50 disabled:cursor-not-allowed transition-colors"
          >
            {submitting ? (
              <>
                <Loader2 className="w-5 h-5 animate-spin" />
                Submitting...
              </>
            ) : (
              <>
                <Upload className="w-5 h-5" />
                Submit Job
              </>
            )}
          </button>
        </div>
      </form>
    </div>
  )
}
