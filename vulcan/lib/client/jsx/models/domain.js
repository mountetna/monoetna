export default class Domain {
  constructor(data) {
    this.data = data;
  }

  get asArray() {
    return [this.padInt(Math.min), this.padInt(Math.max)];
  }

  padInt(func) {
    const margin = 1.2;
    return Math.round(margin * func(...this.data));
  }
}
