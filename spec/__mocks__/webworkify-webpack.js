import {Debouncer} from "../../utils/debouncer";

export default function mockWebworkify(script) {
  const outHandlers = [];
  const inHandlers = [];
  const outQueue = [];
  const inQueue = [];
  const debouncer = new Debouncer();
  let terminated = false;

  const setup = require(script).default

  setTimeout(() => setup({
    postMessage(data) {
      outQueue.push(data);
      debouncer.ready(() => {
        if (terminated) return;
        const queue = outQueue.slice();
        outQueue.length = 0;
        queue.forEach(data => outHandlers.forEach(handler => handler({ data })));
      })
    },
    addEventListener(type, handler) {
      if (type === 'message') inHandlers.push(handler);
    },
  }), 0);

  return {
    postMessage(data) {
      inQueue.push(data);
      debouncer.ready(() => {
        if (terminated) return;
        const queue = inQueue.slice();
        inQueue.length = 0;
        queue.forEach(data => inHandlers.forEach(handler => handler({ data })));
      })
    },
    terminate() {
      terminated = true;
    },
    addEventListener(type, handler) {
      if (type === 'message') outHandlers.push(handler);
    },
    removeEventListener(type, handler) {
      if (type === 'message' && outHandlers.includes(handler)) outHandlers.splice(outHandlers.indexOf(handler), 1);
    },
  }
}