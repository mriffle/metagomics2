import { useCallback, useEffect, useRef, useState } from 'react'
import { X } from 'lucide-react'

export interface AutocompleteOption {
  value: string
  label: string
}

interface AutocompleteProps {
  options: AutocompleteOption[]
  value: string
  onChange: (value: string) => void
  placeholder?: string
  maxResults?: number
}

export default function Autocomplete({
  options,
  value,
  onChange,
  placeholder = 'Search...',
  maxResults = 50,
}: AutocompleteProps) {
  const [inputText, setInputText] = useState('')
  const [isOpen, setIsOpen] = useState(false)
  const [highlightIndex, setHighlightIndex] = useState(-1)
  const containerRef = useRef<HTMLDivElement>(null)
  const inputRef = useRef<HTMLInputElement>(null)
  const listRef = useRef<HTMLUListElement>(null)

  // Sync input text when external value changes (e.g. cleared)
  useEffect(() => {
    if (value) {
      const opt = options.find(o => o.value === value)
      setInputText(opt ? opt.label : value)
    } else {
      setInputText('')
    }
  }, [value, options])

  const filtered = inputText && !value
    ? options
        .filter(o => o.label.toLowerCase().includes(inputText.toLowerCase()))
        .slice(0, maxResults)
    : []

  const handleInputChange = (text: string) => {
    setInputText(text)
    setHighlightIndex(-1)
    if (value) {
      // User is editing after a selection — clear the selection
      onChange('')
    }
    setIsOpen(text.length > 0)
  }

  const handleSelect = useCallback((optValue: string) => {
    onChange(optValue)
    setIsOpen(false)
    setHighlightIndex(-1)
    inputRef.current?.blur()
  }, [onChange])

  const handleClear = () => {
    setInputText('')
    onChange('')
    setIsOpen(false)
    setHighlightIndex(-1)
    inputRef.current?.focus()
  }

  const handleKeyDown = (e: React.KeyboardEvent) => {
    if (!isOpen || filtered.length === 0) return

    if (e.key === 'ArrowDown') {
      e.preventDefault()
      setHighlightIndex(prev => Math.min(prev + 1, filtered.length - 1))
    } else if (e.key === 'ArrowUp') {
      e.preventDefault()
      setHighlightIndex(prev => Math.max(prev - 1, 0))
    } else if (e.key === 'Enter') {
      e.preventDefault()
      if (highlightIndex >= 0 && highlightIndex < filtered.length) {
        handleSelect(filtered[highlightIndex].value)
      }
    } else if (e.key === 'Escape') {
      setIsOpen(false)
      setHighlightIndex(-1)
    }
  }

  // Scroll highlighted item into view
  useEffect(() => {
    if (highlightIndex >= 0 && listRef.current) {
      const item = listRef.current.children[highlightIndex] as HTMLElement
      item?.scrollIntoView({ block: 'nearest' })
    }
  }, [highlightIndex])

  // Close dropdown on outside click
  useEffect(() => {
    const handleClickOutside = (e: MouseEvent) => {
      if (containerRef.current && !containerRef.current.contains(e.target as Node)) {
        setIsOpen(false)
        setHighlightIndex(-1)
      }
    }
    document.addEventListener('mousedown', handleClickOutside)
    return () => document.removeEventListener('mousedown', handleClickOutside)
  }, [])

  return (
    <div ref={containerRef} className="relative">
      <div className="flex items-center">
        <input
          ref={inputRef}
          type="text"
          value={inputText}
          onChange={(e) => handleInputChange(e.target.value)}
          onFocus={() => { if (inputText && !value) setIsOpen(true) }}
          onKeyDown={handleKeyDown}
          placeholder={placeholder}
          className="w-56 text-xs border border-gray-300 dark:border-gray-600 rounded px-2 py-1 bg-white dark:bg-gray-800 text-gray-900 dark:text-gray-100 placeholder-gray-400 dark:placeholder-gray-500 focus:outline-none focus:ring-1 focus:ring-indigo-500"
        />
        {(inputText || value) && (
          <button
            onClick={handleClear}
            className="ml-1 p-0.5 text-gray-400 dark:text-gray-500 hover:text-gray-600 dark:hover:text-gray-300 transition-colors"
            title="Clear filter"
          >
            <X className="w-3.5 h-3.5" />
          </button>
        )}
      </div>

      {isOpen && filtered.length > 0 && (
        <ul
          ref={listRef}
          className="absolute z-50 mt-1 w-80 max-h-60 overflow-auto bg-white dark:bg-gray-800 border border-gray-200 dark:border-gray-600 rounded-md shadow-lg dark:shadow-indigo-500/10"
        >
          {filtered.map((opt, idx) => (
            <li
              key={opt.value}
              onMouseDown={(e) => { e.preventDefault(); handleSelect(opt.value) }}
              onMouseEnter={() => setHighlightIndex(idx)}
              className={`px-3 py-1.5 text-xs cursor-pointer ${
                idx === highlightIndex
                  ? 'bg-indigo-50 dark:bg-indigo-900/50 text-indigo-700 dark:text-indigo-300'
                  : 'text-gray-700 dark:text-gray-300 hover:bg-gray-50 dark:hover:bg-gray-700'
              }`}
            >
              {opt.label}
            </li>
          ))}
        </ul>
      )}
    </div>
  )
}
