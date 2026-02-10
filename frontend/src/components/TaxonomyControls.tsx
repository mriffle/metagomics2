import { useState } from 'react'
import { Download, Filter } from 'lucide-react'
import type { ChartType } from '../pages/TaxonomyPage'
import { CANONICAL_RANKS_ORDERED } from '../utils/taxonomyParser'
import type { CanonicalRank } from '../utils/taxonomyParser'
import Autocomplete from './Autocomplete'
import type { AutocompleteOption } from './Autocomplete'

interface TaxonomyControlsProps {
  chartType: ChartType
  onChartTypeChange: (t: ChartType) => void
  maxRank: CanonicalRank
  onMaxRankChange: (r: CanonicalRank) => void
  minRatioTotal: number
  onMinRatioTotalChange: (v: number) => void
  filteredNodeCount: number
  totalNodeCount: number
  onExportPng?: () => void
  onExportSvg?: () => void
  goOptions?: AutocompleteOption[]
  selectedGoTerm?: string
  onGoTermChange?: (goId: string) => void
}

const CUTOFF_PRESETS = [
  { label: 'None', value: 0 },
  { label: '0.01%', value: 0.0001 },
  { label: '0.1%', value: 0.001 },
  { label: '1%', value: 0.01 },
  { label: '5%', value: 0.05 },
  { label: '10%', value: 0.1 },
]

const CHART_TYPES: { value: ChartType; label: string }[] = [
  { value: 'sunburst', label: 'Sunburst' },
  { value: 'treemap', label: 'Treemap' },
  { value: 'icicle', label: 'Icicle' },
  { value: 'sankey', label: 'Sankey' },
]

// Selectable ranks (exclude root — it's always included implicitly)
const SELECTABLE_RANKS = CANONICAL_RANKS_ORDERED.filter(r => r !== 'root')

const RANK_LABELS: Record<string, string> = {
  domain: 'Domain',
  kingdom: 'Kingdom',
  phylum: 'Phylum',
  class: 'Class',
  order: 'Order',
  family: 'Family',
  genus: 'Genus',
  species: 'Species',
}

export default function TaxonomyControls({
  chartType,
  onChartTypeChange,
  maxRank,
  onMaxRankChange,
  minRatioTotal,
  onMinRatioTotalChange,
  filteredNodeCount,
  totalNodeCount,
  onExportPng,
  onExportSvg,
  goOptions,
  selectedGoTerm,
  onGoTermChange,
}: TaxonomyControlsProps) {
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
      <div className="flex items-center justify-between flex-wrap gap-2">
      {/* Left: chart type toggle + rank selector */}
      <div className="flex items-center gap-4">
        <div className="flex items-center gap-1">
          {CHART_TYPES.map((ct) => (
            <button
              key={ct.value}
              onClick={() => onChartTypeChange(ct.value)}
              className={`px-3 py-1.5 text-sm rounded-md transition-colors ${
                ct.value === chartType
                  ? 'bg-indigo-600 text-white'
                  : 'bg-gray-100 text-gray-700 hover:bg-gray-200'
              }`}
            >
              {ct.label}
            </button>
          ))}
        </div>

        {/* Rank depth selector */}
        <div className="flex items-center gap-1.5">
          <label htmlFor="max-rank" className="text-sm text-gray-600">Depth:</label>
          <select
            id="max-rank"
            value={maxRank}
            onChange={(e) => onMaxRankChange(e.target.value as CanonicalRank)}
            className="text-sm border border-gray-300 rounded-md px-2 py-1 bg-white focus:outline-none focus:ring-2 focus:ring-indigo-500"
          >
            {SELECTABLE_RANKS.map((r) => (
              <option key={r} value={r}>{RANK_LABELS[r]}</option>
            ))}
          </select>
        </div>
      </div>

      {/* Right: export */}
      <div className="flex items-center gap-3">
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

      {/* GO term filter row */}
      {goOptions && goOptions.length > 0 && onGoTermChange && (
        <div className="flex items-center gap-2 flex-wrap">
          <Filter className="w-3.5 h-3.5 text-gray-500" />
          <span className="text-xs text-gray-600 font-medium">Filter by GO term:</span>
          <Autocomplete
            options={goOptions}
            value={selectedGoTerm || ''}
            onChange={onGoTermChange}
            placeholder="Search GO term..."
          />
          {selectedGoTerm && (
            <span className="text-xs text-indigo-600 font-medium">
              Showing taxonomy for {selectedGoTerm}
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
        <span className="ml-2 text-xs text-gray-500">
          {filteredNodeCount === totalNodeCount
            ? `${filteredNodeCount} nodes`
            : `${filteredNodeCount} / ${totalNodeCount} nodes`}
        </span>
      </div>
    </div>
  )
}
