// Series for the data object, not the plot definition
import Gradient from 'javascript-color-gradient';
import Vector from 'etna-js/plots/models/vector';
import {autoColors} from 'etna-js/utils/colors';

export default class Series {
  public name: string | null = null;
  public type: string | null = null;
  public xValues: any[] = [];
  public yValues: any[] = [];
  public labels: string[] = [];
  public secondaryData: {rows: number[][], row_names: string[]} | null = null;
  public colorBy: string | null = null;

  constructor() {}

  get asObj() {
    return {
      name: this.name,
      series_type: this.type,
      variables: {
        x: this.vectorize(this.xValues),
        y: this.vectorize(this.yValues),
        label: this.vectorize(this.labels),
        nodeColor: this.colorBy ? this.colors() : null
      }
    };
  }

  setName(name: string | null) {
    this.name = name;
  }

  setType(type: string | null) {
    this.type = type;
  }

  setXValues(xValues: any[]) {
    this.xValues = xValues;
  }

  setYValues(yValues: any[]) {
    this.yValues = yValues;
  }

  setLabels(labels: string[]) {
    this.labels = labels;
  }

  setSecondaryData(dataMatrix: Series['secondaryData']) {
    // Assumes rows are selectable options and
    //   columns correspond to the x and y values.
    this.secondaryData = dataMatrix;
  }

  vectorize(values: any[]) {
    return new Vector(
      values.map((val) => ({
        label: null,
        value: val
      }))
    );
  }

  setColorBy(val: string | null) {
    this.colorBy = val ? val : null;
  }

  colors() {
    // Sets the colors based on the secondaryData
    if (
      !this.secondaryData ||
      !this.colorBy ||
      -1 === this.secondaryData.row_names.indexOf(this.colorBy)
    )
      return null;

    let [startColor, endColor] = autoColors(2);

    const colorGradient = new Gradient();

    // Divide up the range by 100, and then set
    //   the colors based on percentage.
    const rowIndex = this.secondaryData.row_names.indexOf(this.colorBy);
    const slice = this.secondaryData.rows[rowIndex];
    if (!slice) return null;

    const minValue = Math.min(...slice);
    const maxValue = Math.max(...slice);

    colorGradient.setMidpoint(100);
    colorGradient.setGradient(startColor || '', endColor || '');

    return slice.map((val: number) => {
      return colorGradient.getColor(
        Math.round((100 * (val - minValue)) / (maxValue - minValue)) + 1
      );
    });
  }
}
