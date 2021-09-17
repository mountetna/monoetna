import {act, ReactTestInstance} from "react-test-renderer";

export function findNode(node: ReactTestInstance, predicate: (e: ReactTestInstance) => boolean, idx: number) {
  const nodes = node.findAll(predicate);
  if (idx < 0) idx = nodes.length + idx;
  return nodes[idx];
}

export async function clickNode(node: ReactTestInstance, predicate: (e: ReactTestInstance) => boolean, idx: number) {
  await act(async () => {
    findNode(node, predicate, idx).props.onClick();
  })
}

export function includesClassNamePredicate(className: string) {
  return (node: ReactTestInstance) => node.props.className?.match(new RegExp(`(?:\\s+|^)${className}(?:\\s+|$)`))
}

export function matchesTypePredicate(element: string) {
  return (node: ReactTestInstance) => node.type === element;
}

export function matchesTextPredicate(t: string) {
  return (node: ReactTestInstance) => text(node) === t;
}

export function findAllByClassName(element: ReactTestInstance, className: string) {
  return element.findAll(includesClassNamePredicate(className));
}

export function text(element: ReactTestInstance) {
  let result = "";
  element.children?.forEach(node => {
    if (typeof node === "string") result += node;
    else result += text(node);
  })

  return result;
}