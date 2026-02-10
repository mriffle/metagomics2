import { useState } from 'react'
import { Download, Filter } from 'lucide-react'
import type { MetricKey } from '../pages/GoDagPage'
import Autocomplete from './Autocomplete'
import type { AutocompleteOption } from './Autocomplete'

interface GoDagControlsProps {
  namespaces: string[]
  namespaceLabels: Record<string, string>
  selectedNamespace: string
  onNamespaceChange: (ns: string) => void
  selectedMetric: MetricKey
  onMetricChange: (m: MetricKey) => void
  minRatioTotal: number
  onMinRatioTotalChange: (v: number) => void
  filteredNodeCount: number
  totalNodeCount: number
  onExportPng?: () => void
  onExportSvg?: () => void
  taxonOptions?: AutocompleteOption[]
  selectedTaxon?: string
  onTaxonChange?: (taxId: string) => void
  filterLabel?: string
}

const CUTOFF_PRESETS = [
  { label: 'None', value: 0 },
  { label: '0.01%', value: 0.0001 },
  { label: '0.1%', value: 0.001 },
  { label: '1%', value: 0.01 },
  { label: '5%', value: 0.05 },
  { label: '10%', value: 0.1 },
]

const BASE_METRIC_OPTIONS: { value: MetricKey; label: string }[] = [
  { value: 'quantity', label: 'Quantity' },
  { value: 'ratioTotal', label: 'Ratio (Total)' },
  { value: 'nPeptides', label: '# Peptides' },
]

const FILTER_METRIC_OPTIONS: { value: MetricKey; label: string }[] = [
  { value: 'fractionOfTaxon', label: 'Fraction of Taxon' },
  { value: 'fractionOfGo', label: 'Fraction of GO' },
]

export default function GoDagControls({
  namespaces,
  namespaceLabels,
  selectedNamespace,
  onNamespaceChange,
  selectedMetric,
  onMetricChange,
  minRatioTotal,
  onMinRatioTotalChange,
  filteredNodeCount,
  totalNodeCount,
  onExportPng,
  onExportSvg,
  taxonOptions,
  selectedTaxon,
  onTaxonChange,
  filterLabel,
}: GoDagControlsProps) {
  const [customCutoff, setCustomCutoff] = useState('')

  const handleCustomCutoff = () => {
    const val = parseFloat(customCutoff)
    if (!isNaN(val) && val >= 0 && val <= 100) {
      onMinRatioTotalChange(val / 100)
      setCustomCutoff('')
    }
  }

  return (
    <div className="flex flex-col gap-2 mb-3">
      {/* Top row: namespace tabs + right-side controls */}
      <div className="flex items-center justify-between flex-wrap gap-2">
        {/* Namespace tabs */}
        <div className="flex items-center gap-1">
          {namespaces.map((ns) => (
            <button
              key={ns}
              onClick={() => onNamespaceChange(ns)}
              className={`px-3 py-1.5 text-sm rounded-md transition-colors ${
                ns === selectedNamespace
                  ? 'bg-indigo-600 text-white'
                  : 'bg-gray-100 text-gray-700 hover:bg-gray-200'
              }`}
            >
              {namespaceLabels[ns] || ns}
            </button>
          ))}
          <span className="ml-2 text-xs text-gray-500">
            {filteredNodeCount === totalNodeCount
              ? `${filteredNodeCount} terms`
              : `${filteredNodeCount} / ${totalNodeCount} terms`}
          </span>
        </div>

        {/* Right side controls */}
        <div className="flex items-center gap-3">
          {/* Metric selector */}
          <div className="flex items-center gap-2">
            <label className="text-xs text-gray-600 font-medium">Color by:</label>
            <select
              value={selectedMetric}
              onChange={(e) => onMetricChange(e.target.value as MetricKey)}
              className="text-sm border border-gray-300 rounded-md px-2 py-1 bg-white focus:outline-none focus:ring-1 focus:ring-indigo-500"
            >
              {(filterLabel ? [...BASE_METRIC_OPTIONS, ...FILTER_METRIC_OPTIONS] : BASE_METRIC_OPTIONS).map((opt) => (
                <option key={opt.value} value={opt.value}>
                  {opt.label}
                </option>
              ))}
            </select>
          </div>

          {/* Color legend */}
          <div className="flex items-center gap-1">
            <span className="text-xs text-gray-500">Low</span>
            <div
              className="w-20 h-3 rounded"
              style={{
                background: 'linear-gradient(to right, #eef2ff, #4338ca)',
              }}
            />
            <span className="text-xs text-gray-500">High</span>
          </div>

          {/* Export */}
          {onExportPng && (
            <button
              onClick={onExportPng}
              className="inline-flex items-center gap-1 px-2 py-1 text-xs bg-gray-100 hover:bg-gray-200 rounded transition-colors"
            >
              <Download className="w-3 h-3" />
              PNG
            </button>
          )}
          {onExportSvg && (
            <button
              onClick={onExportSvg}
              className="inline-flex items-center gap-1 px-2 py-1 text-xs bg-gray-100 hover:bg-gray-200 rounded transition-colors"
            >
              <Download className="w-3 h-3" />
              SVG
            </button>
          )}
        </div>
      </div>

      {/* Taxonomy filter row */}
      {taxonOptions && taxonOptions.length > 0 && onTaxonChange && (
        <div className="flex items-center gap-2 flex-wrap">
          <Filter className="w-3.5 h-3.5 text-gray-500" />
          <span className="text-xs text-gray-600 font-medium">Filter by taxon:</span>
          <Autocomplete
            options={taxonOptions}
            value={selectedTaxon || ''}
            onChange={onTaxonChange}
            placeholder="Search taxonomy..."
          />
          {selectedTaxon && (
            <span className="text-xs text-indigo-600 font-medium">
              Showing GO terms for taxon {selectedTaxon}
            </span>
          )}
        </div>
      )}

      {/* Bottom row: abundance cutoff filter */}
      <div className="flex items-center gap-2 flex-wrap">
        <Filter className="w-3.5 h-3.5 text-gray-500" />
        <span className="text-xs text-gray-600 font-medium">Min ratio (total):</span>
        {CUTOFF_PRESETS.map((preset) => (
          <button
            key={preset.value}
            onClick={() => onMinRatioTotalChange(preset.value)}
            className={`px-2 py-0.5 text-xs rounded transition-colors ${
              minRatioTotal === preset.value
                ? 'bg-indigo-600 text-white'
                : 'bg-gray-100 text-gray-700 hover:bg-gray-200'
            }`}
          >
            {preset.label}
          </button>
        ))}
        <div className="flex items-center gap-1">
          <input
            type="text"
            value={customCutoff}
            onChange={(e) => setCustomCutoff(e.target.value)}
            onKeyDown={(e) => e.key === 'Enter' && handleCustomCutoff()}
            placeholder="%"
            className="w-14 text-xs border border-gray-300 rounded px-1.5 py-0.5 focus:outline-none focus:ring-1 focus:ring-indigo-500"
          />
          <button
            onClick={handleCustomCutoff}
            className="px-1.5 py-0.5 text-xs bg-gray-100 hover:bg-gray-200 rounded transition-colors"
          >
            Set
          </button>
        </div>
        {minRatioTotal > 0 && !CUTOFF_PRESETS.some(p => p.value === minRatioTotal) && (
          <span className="text-xs text-indigo-600 font-medium">
            {(minRatioTotal * 100).toFixed(4).replace(/0+$/, '').replace(/\.$/, '')}%
          </span>
        )}
      </div>
    </div>
  )
}
