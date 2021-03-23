export default class Domain {
  constructor(public data: number[]) {
  }

  get asArray() {
    return [this.padInt(Math.min), this.padInt(Math.max)];
  }

  padInt(func: (...t: number[]) => number) {
    const margin = 1.2;
    return Math.round(margin * func(...this.data));
  }
};
