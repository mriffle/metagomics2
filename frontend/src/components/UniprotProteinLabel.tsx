import { useEffect, useLayoutEffect, useRef, useState } from 'react'
import { ExternalLink, Loader2, ShieldCheck, ShieldAlert } from 'lucide-react'
import { uniprotUrl, fetchUniprotInfo, extractUniprotAccession } from '../utils/uniprot'
import type { UniprotInfo } from '../utils/uniprot'

interface TooltipContentProps {
  info: UniprotInfo
  accession: string
}

function TooltipContent({ info, accession }: TooltipContentProps) {
  return (
    <div className="text-xs space-y-1">
      <div className="flex items-center gap-1.5 font-semibold text-gray-900">
        {info.reviewed
          ? <ShieldCheck className="w-3.5 h-3.5 text-amber-500 flex-shrink-0" />
          : <ShieldAlert className="w-3.5 h-3.5 text-gray-400 flex-shrink-0" />}
        <span>{info.reviewed ? 'Swiss-Prot (reviewed)' : 'TrEMBL (unreviewed)'}</span>
      </div>
      <div className="text-gray-800 leading-snug">{info.fullName}</div>
      {info.gene && (
        <div className="text-gray-500">Gene: <span className="font-mono text-gray-700">{info.gene}</span></div>
      )}
      {info.organism && (
        <div className="text-gray-500 italic">{info.organism}</div>
      )}
      <div className="pt-0.5 text-indigo-600 font-mono">{accession}</div>
    </div>
  )
}

interface UniprotProteinLabelProps {
  rawId: string
}

export default function UniprotProteinLabel({ rawId }: UniprotProteinLabelProps) {
  const accession = extractUniprotAccession(rawId)!
  const [info, setInfo] = useState<UniprotInfo | null | 'loading' | 'error'>('loading')
  const [showTooltip, setShowTooltip] = useState(false)
  const anchorRef = useRef<HTMLAnchorElement>(null)
  const tooltipRef = useRef<HTMLDivElement>(null)
  const [tooltipPos, setTooltipPos] = useState<{ left: number; top: number }>({ left: 0, top: 0 })
  const fetchedRef = useRef(false)

  useEffect(() => {
    if (fetchedRef.current) return
    fetchedRef.current = true
    fetchUniprotInfo(accession).then(result => {
      setInfo(result ?? 'error')
    })
  }, [accession])

  useLayoutEffect(() => {
    if (!showTooltip || !anchorRef.current || !tooltipRef.current) return
    const aRect = anchorRef.current.getBoundingClientRect()
    const tRect = tooltipRef.current.getBoundingClientRect()

    let left = aRect.left
    let top = aRect.bottom + 4

    if (left + tRect.width > window.innerWidth - 8) {
      left = window.innerWidth - tRect.width - 8
    }
    if (left < 8) left = 8

    if (top + tRect.height > window.innerHeight - 8) {
      top = aRect.top - tRect.height - 4
    }

    setTooltipPos({ left, top })
  }, [showTooltip, info])

  const hasInfo = info !== 'loading' && info !== 'error' && info !== null

  return (
    <span className="relative inline-flex items-center gap-1">
      <a
        ref={anchorRef}
        href={uniprotUrl(accession)}
        target="_blank"
        rel="noopener noreferrer"
        className="font-mono text-xs text-emerald-700 hover:text-emerald-900 hover:underline break-all"
        onMouseEnter={() => setShowTooltip(true)}
        onMouseLeave={() => setShowTooltip(false)}
        onClick={e => e.stopPropagation()}
      >
        {accession}
      </a>

      {info === 'loading' && (
        <Loader2 className="w-2.5 h-2.5 text-gray-300 animate-spin flex-shrink-0" />
      )}

      {hasInfo && (
        <ExternalLink className="w-2.5 h-2.5 text-gray-300 flex-shrink-0" />
      )}

      {showTooltip && hasInfo && (
        <div
          ref={tooltipRef}
          className="fixed z-50 w-64 bg-white border border-gray-200 rounded-lg shadow-lg p-3 pointer-events-none"
          style={{ left: tooltipPos.left, top: tooltipPos.top }}
        >
          <TooltipContent info={info as UniprotInfo} accession={accession} />
        </div>
      )}
    </span>
  )
}
