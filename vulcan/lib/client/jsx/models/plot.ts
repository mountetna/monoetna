import {WorkflowStep} from "../api_types";
import Series from "./series";

export default class Plot {
  public configObj: any = {};
  public dataObj: any = {};
  public series: Series[] = [];


  constructor(public step: WorkflowStep, public consignment: any, public parentWidth: number) {
  }

  get config() {
    return this.configObj;
  }

  get data() {
    return this.dataObj;
  }

  get validSeries() {
    // Because we may have null series in the List,
    //   we'll need to filter those out and only
    //   use the "real" series.
    return this.series.filter((s) => null !== s);
  }

  get hasColorableSeries() {
    // Check if any valid series has secondaryData attached
    return this.validSeries.some((s) => s.secondaryData);
  }

  get minsMaxes() {
    // Calculate this across all series
    let current = {
      x_min: Math.min(), // Counterintuitively, returns Infinity
      x_max: Math.max(), // Counterintuitively, returns -Infinity
      y_min: Math.min(),
      y_max: Math.max()
    };
    this.validSeries.forEach((series) => {
      current.x_min = Math.min(...[...series.xValues, current.x_min]);
      current.x_max = Math.max(...[...series.xValues, current.x_max]);
      current.y_min = Math.min(...[...series.yValues, current.y_min]);
      current.y_max = Math.max(...[...series.yValues, current.y_max]);
    });
    return current;
  }
}
