import { describe, it, expect, beforeEach } from 'vitest'
import { render, screen, fireEvent } from '@testing-library/react'
import { ThemeProvider, useTheme } from '../ThemeContext'

// Helper component that exposes theme context values for testing
function ThemeDisplay() {
  const { theme, toggleTheme } = useTheme()
  return (
    <div>
      <span data-testid="theme-value">{theme}</span>
      <button onClick={toggleTheme}>Toggle</button>
    </div>
  )
}

describe('ThemeContext', () => {
  beforeEach(() => {
    // Clear localStorage and remove dark class before each test
    localStorage.clear()
    document.documentElement.classList.remove('dark')
  })

  it('defaults to light theme when no localStorage value', () => {
    render(
      <ThemeProvider>
        <ThemeDisplay />
      </ThemeProvider>
    )
    expect(screen.getByTestId('theme-value').textContent).toBe('light')
    expect(document.documentElement.classList.contains('dark')).toBe(false)
  })

  it('reads initial theme from localStorage', () => {
    localStorage.setItem('metagomics-theme', 'dark')
    render(
      <ThemeProvider>
        <ThemeDisplay />
      </ThemeProvider>
    )
    expect(screen.getByTestId('theme-value').textContent).toBe('dark')
    expect(document.documentElement.classList.contains('dark')).toBe(true)
  })

  it('toggles from light to dark', () => {
    render(
      <ThemeProvider>
        <ThemeDisplay />
      </ThemeProvider>
    )
    expect(screen.getByTestId('theme-value').textContent).toBe('light')

    fireEvent.click(screen.getByText('Toggle'))

    expect(screen.getByTestId('theme-value').textContent).toBe('dark')
    expect(document.documentElement.classList.contains('dark')).toBe(true)
    expect(localStorage.getItem('metagomics-theme')).toBe('dark')
  })

  it('toggles from dark to light', () => {
    localStorage.setItem('metagomics-theme', 'dark')
    render(
      <ThemeProvider>
        <ThemeDisplay />
      </ThemeProvider>
    )
    expect(screen.getByTestId('theme-value').textContent).toBe('dark')

    fireEvent.click(screen.getByText('Toggle'))

    expect(screen.getByTestId('theme-value').textContent).toBe('light')
    expect(document.documentElement.classList.contains('dark')).toBe(false)
    expect(localStorage.getItem('metagomics-theme')).toBe('light')
  })

  it('persists theme to localStorage on toggle', () => {
    render(
      <ThemeProvider>
        <ThemeDisplay />
      </ThemeProvider>
    )

    // Toggle to dark
    fireEvent.click(screen.getByText('Toggle'))
    expect(localStorage.getItem('metagomics-theme')).toBe('dark')

    // Toggle back to light
    fireEvent.click(screen.getByText('Toggle'))
    expect(localStorage.getItem('metagomics-theme')).toBe('light')
  })

  it('ignores invalid localStorage values and defaults to light', () => {
    localStorage.setItem('metagomics-theme', 'invalid')
    render(
      <ThemeProvider>
        <ThemeDisplay />
      </ThemeProvider>
    )
    expect(screen.getByTestId('theme-value').textContent).toBe('light')
  })

  it('adds dark class to documentElement in dark mode', () => {
    render(
      <ThemeProvider>
        <ThemeDisplay />
      </ThemeProvider>
    )

    fireEvent.click(screen.getByText('Toggle'))
    expect(document.documentElement.classList.contains('dark')).toBe(true)

    fireEvent.click(screen.getByText('Toggle'))
    expect(document.documentElement.classList.contains('dark')).toBe(false)
  })
})
