// Series for the data object, not the plot definition
import Gradient from 'javascript-color-gradient';

import Vector from 'etna-js/plots/models/vector';
import {autoColors} from 'etna-js/utils/colors';

export default class Series {
  constructor() {
    this.name = null;
    this.type = null;
    this.xValues = null;
    this.yValues = null;
    this.labels = null;
    this.secondaryData = null;
    this.colorBy = null;
  }

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

  setName(name) {
    this.name = name;
  }

  setType(type) {
    this.type = type;
  }

  setXValues(xValues) {
    this.xValues = xValues;
  }

  setYValues(yValues) {
    this.yValues = yValues;
  }

  setLabels(labels) {
    this.labels = labels;
  }

  setSecondaryData(dataMatrix) {
    // Assumes rows are selectable options and
    //   columns correspond to the x and y values.
    this.secondaryData = dataMatrix;
  }

  vectorize(values) {
    return new Vector(
      values.map((val) => ({
        label: null,
        value: val
      }))
    );
  }

  setColorBy(val) {
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

    let colorOptions = autoColors(2);
    let startColor = colorOptions[0];
    let endColor = colorOptions[1];

    const colorGradient = new Gradient();

    // Divide up the range by 100, and then set
    //   the colors based on percentage.
    const rowIndex = this.secondaryData.row_names.indexOf(this.colorBy);
    const slice = this.secondaryData.rows[rowIndex];
    const minValue = Math.min(...slice);
    const maxValue = Math.max(...slice);

    colorGradient.setMidpoint(100);
    colorGradient.setGradient(startColor, endColor);

    return slice.map((val) => {
      return colorGradient.getColor(
        Math.round((100 * (val - minValue)) / (maxValue - minValue)) + 1
      );
    });
  }
}
