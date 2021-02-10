export default class Plot {
  constructor(step, consignment, dimensions) {
    this.step = step;
    this.consignment = consignment;
    this.dimensions = dimensions;

    this.plotObj = {};
    this.dataObj = {};
  }

  get plot() {
    return this.plotObj;
  }

  get data() {
    return this.dataObj;
  }
}
