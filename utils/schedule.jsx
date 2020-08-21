// A schedule is simply a continually resizing group of pending promises that can be
// awaited as a group.
export class Schedule {
  constructor() {
    this.work = [];
  }

  addWork(work) {
    work = work.finally(() => this.work = this.work.filter(v => v !== work));
    this.work.push(work);
    return work;
  }

  curPending() {
    return Promise.all(this.work);
  }

  allPending() {
    return this.curPending().then(v => {
      if (this.work.length > 0) return this.curPending().then(vv => v.concat(vv));
      return v;
    })
  }
}