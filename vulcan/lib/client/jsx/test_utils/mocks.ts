export function createFakeStorage(): Storage {
  const storage: { [k: string]: string } = {};

  return {
    setItem(key: string, value: string) {
      storage[key] = value + "";
    },
    get length() {
      return Object.keys(storage).length;
    },
    clear() {
      Object.keys(storage).forEach(k => delete storage[k]);
    },
    getItem(key: string): string | null {
      if (key in storage) return storage[key];
      return null;
    },
    removeItem(key: string) {
      delete storage[key];
    },
    key(index: number): string | null {
      return Object.keys(storage)[index];
    }
  }
}