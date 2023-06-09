
export const diff = (list1: string[], list2: string[]) => {
  const l2 = new Set(list2);
  return list1.filter((x) => !l2.has(x));
};

