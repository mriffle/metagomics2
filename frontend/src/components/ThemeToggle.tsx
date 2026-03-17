import { Sun, Moon } from 'lucide-react'
import { useTheme } from '../ThemeContext'

export default function ThemeToggle() {
  const { theme, toggleTheme } = useTheme()
  const isDark = theme === 'dark'

  return (
    <button
      onClick={toggleTheme}
      aria-label={isDark ? 'Switch to light mode' : 'Switch to dark mode'}
      className="relative inline-flex h-7 w-14 items-center rounded-full transition-colors duration-300 focus:outline-none focus:ring-2 focus:ring-indigo-500 focus:ring-offset-2 dark:focus:ring-offset-gray-900 bg-gray-200 dark:bg-indigo-600"
    >
      {/* Sliding knob */}
      <span
        className={`inline-flex h-5 w-5 items-center justify-center rounded-full bg-white shadow-md transition-transform duration-300 ${
          isDark ? 'translate-x-8' : 'translate-x-1'
        }`}
      >
        {isDark ? (
          <Moon className="h-3 w-3 text-indigo-600" />
        ) : (
          <Sun className="h-3 w-3 text-amber-500" />
        )}
      </span>

      {/* Background icons */}
      <Sun className={`absolute left-1.5 h-3.5 w-3.5 transition-opacity duration-300 ${isDark ? 'opacity-30 text-gray-400' : 'opacity-0'}`} />
      <Moon className={`absolute right-1.5 h-3.5 w-3.5 transition-opacity duration-300 ${isDark ? 'opacity-0' : 'opacity-30 text-gray-500'}`} />
    </button>
  )
}
