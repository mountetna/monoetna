export default class Plot {
  constructor(step, consignment, parentWidth) {
    this.step = step;
    this.consignment = consignment;
    this.parentWidth = parentWidth;

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
