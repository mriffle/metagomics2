import { Download } from 'lucide-react'
import type { ChartType } from '../pages/TaxonomyPage'
import { CANONICAL_RANKS_ORDERED } from '../utils/taxonomyParser'
import type { CanonicalRank } from '../utils/taxonomyParser'

interface TaxonomyControlsProps {
  chartType: ChartType
  onChartTypeChange: (t: ChartType) => void
  maxRank: CanonicalRank
  onMaxRankChange: (r: CanonicalRank) => void
  onExportPng?: () => void
  onExportSvg?: () => void
}

const CHART_TYPES: { value: ChartType; label: string }[] = [
  { value: 'sunburst', label: 'Sunburst' },
  { value: 'treemap', label: 'Treemap' },
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
  onExportPng,
  onExportSvg,
}: TaxonomyControlsProps) {
  return (
    <div className="flex items-center justify-between flex-wrap gap-2 mb-3">
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
  )
}
