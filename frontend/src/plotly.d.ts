declare module 'plotly.js-dist-min' {
  import Plotly from 'plotly.js'
  export default Plotly
  export * from 'plotly.js'
}

declare module 'react-plotly.js/factory' {
  import { Component } from 'react'
  import Plotly from 'plotly.js'

  interface PlotParams {
    data: Plotly.Data[]
    layout?: Partial<Plotly.Layout>
    config?: Partial<Plotly.Config>
    frames?: Plotly.Frame[]
    style?: React.CSSProperties
    className?: string
    useResizeHandler?: boolean
    onHover?: (event: any) => void
    onUnhover?: (event: any) => void
    onClick?: (event: any) => void
    onSelected?: (event: any) => void
    onRelayout?: (event: any) => void
    onUpdate?: (figure: any, graphDiv: HTMLElement) => void
    onInitialized?: (figure: any, graphDiv: HTMLElement) => void
    onPurge?: (figure: any, graphDiv: HTMLElement) => void
    revision?: number
    divId?: string
  }

  function createPlotlyComponent(plotly: typeof Plotly): React.ComponentType<PlotParams>
  export default createPlotlyComponent
}
