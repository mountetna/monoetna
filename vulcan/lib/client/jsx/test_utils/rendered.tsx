import {ReactTestInstance} from "react-test-renderer";

export function includesClassNamePredicate(className: string) {
  return (node: ReactTestInstance) => node.props.className?.match(new RegExp(`(?:\\s+|^)${className}(?:\\s+|$)`))
}

export function matchesTypePredicate(element: string) {
  return (node: ReactTestInstance) => node.type === element;
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