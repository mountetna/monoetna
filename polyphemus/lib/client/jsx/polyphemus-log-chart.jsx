import React, {useState, useEffect, useCallback, useContext} from 'react';
import Axis from 'etna-js/plots/components/axis';
import Legend from 'etna-js/plots/components/legend';
import * as d3 from 'd3';
import Box from '@material-ui/core/Box';
import { APP_COLORS } from './polyphemus-logs';

const scale = (domain, range) => {
  let s;

  if (domain[0] instanceof Date) {
    s = d3.scaleTime();
    s.type = 'time';
    console.log({s});
  } else if (typeof domain[0] === 'string') {
    s = d3.scaleBand().padding(0.2);
    s.type = 'band';
  } else {
    s = d3.scaleLinear();
    s.type = 'linear';
  }

  s = s.domain(domain).range(range);
  if (s.nice) s = s.nice();

  return s;
}

const PlotCanvas = ({
    labels,
    xdomain,
    ydomain,
    xlabel,
    ylabel,
    children
  }) => {
    console.log({xdomain,ydomain});
  const margin = { left: 40, right: 40, top: 40, bottom: 40 };

  const [stageSize, setStageSize] = useState({width: 100, height: 100});
  const stage = React.useRef(null);

  const resize = () => {
    if (stage.current) setStageSize(stage.current.getBoundingClientRect());
  }

  useEffect( () => {
    resize();
  }, [stage] );

  useEffect(() => {
    window.addEventListener('resize', resize);
    return () => {
      window.removeEventListener('resize', resize);
    };
  }, []);

  let xScale = scale(xdomain, [
    margin.left,
    stageSize.width - margin.right
  ]);

  let yScale = scale(
    ydomain,
    [
      stageSize.height - margin.bottom,
      margin.top
    ]
  );

  let xsize = stageSize.width - margin.left - margin.right;
  let ysize = stageSize.height - margin.top - margin.bottom;

  return (
    <Box
      sx={
        {
          position: 'relative',
          width: '100%', height: '100%',
          '& .tick': {
            color: 'skyblue',
            '& text': {
              color: 'black',
            }
          },
          '& .line': {
            strokeWidth: 5,
            strokeLinecap: 'round',
            strokeLineround: 'round'
          }
        }
      }
      ref={stage}>
      <Legend width={stageSize.width} height={margin.top} labels={labels} />
      <svg width={stageSize.width} height={stageSize.height}>
        <rect
          fill='aliceblue'
          x={ margin.right }
          y={ margin.top }
          width={ stageSize.width - margin.right - margin.left }
          height={ stageSize.height - margin.top - margin.bottom } />
        {xlabel && (
          <text
            x={xsize / 2 + parseInt(margin.left)}
            y={parseInt(margin.top) + ysize + 50}
            fontSize='12px'
            textAnchor='middle'
          >
            {xlabel}
          </text>
        )}
        <Axis
          orient='Bottom'
          scale={xScale}
          translate={`translate(0, ${stageSize.height - margin.bottom})`}
          tickSize={ysize}
        />
        {ylabel && (
          <text
            x={-ysize / 2 - parseInt(margin.top)}
            y={margin.left - 50}
            fontSize='12px'
            transform='rotate(-90)'
            textAnchor='middle'
          >
            {ylabel}
          </text>
        )}
        <Axis
          orient='Left'
          scale={yScale}
          translate={`translate(${margin.left}, 0)`}
          tickSize={xsize}
        />
        {
          React.Children.map(children, child =>
            React.cloneElement(child, { xScale, yScale }))
        }
      </svg>
    </Box>
  );
}

const Line = ({ name, series, xScale, yScale, color }) => {
  const linePath = d3
    .line()
    .x(d => d.x)
    .y(d => d.y);

  const points = series.map((p, i) => ({ x: xScale(p.x), y: yScale(p.y) }));
  const path = linePath(points);

  return <path 
    className="line"
    d={ linePath(points) }
    fill="none"
    stroke={ color }
  />;
}

const HistogramChart = ({logs}) => {
  const fromDate = logs.reduce((min, d) => (min > d.created_at ? d.created_at : min), '3000-01-01');
  const toDate = logs.reduce((max, d) => (max < d.created_at ? d.created_at : max), '1970-01-01');

  let groups = Object.groupBy(logs, log => log.application);
  groups.all = logs

  const xdomain = [new Date(fromDate), new Date(toDate)];

  const xScale = scale(xdomain, [ 0, 1 ]);
  const hist = d3.histogram().value(d => d).domain(xdomain).thresholds(xScale.ticks(d3.timeMonth));

  groups = Object.fromEntries( Object.keys(groups).map(
    group => [ group, {
      data: groups[group],
      bins: hist(groups[group].map( l => new Date(l.created_at)))
    } ]
  ) );

  const maxBinSize = Math.max( ...Object.keys(groups).map( group => Math.max( ...groups[group].bins.map( bin => bin.length ) ) ) );

  console.log({fromDate,toDate, logs, groups, maxBinSize});

  return <PlotCanvas
    xdomain={xdomain}
    ydomain={[-0.5,maxBinSize+1]}
    labels={Object.keys(groups).map(name => ({name, color: APP_COLORS[name] || 'pink' }))}
  >
  {
    Object.keys(groups).map( group => <Line
      key={group}
      name={group}
      color={ APP_COLORS[group] || 'pink' }
      series={
        groups[group].bins.map(
          bin => ({ x: bin.x0, y: bin.length })
        )
      }
    /> )
  }
  </PlotCanvas>;
}

export default HistogramChart;
