// Series for the data object, not the plot definition
import Vector from 'etna-js/plots/models/vector';

export default class Series {
  constructor() {
    this.name = null;
    this.type = null;
    this.xValues = null;
    this.yValues = null;
    this.labels = null;
  }

  get asObj() {
    return {
      name: this.name,
      series_type: this.type,
      variables: {
        x: this.vectorize(this.xValues),
        y: this.vectorize(this.yValues),
        label: this.vectorize(this.labels)
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

  vectorize(values) {
    return new Vector(
      values.map((val) => ({
        label: null,
        value: val
      }))
    );
  }
}
