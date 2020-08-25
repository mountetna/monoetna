// Mock service worker registration that never completes install.
export default function registerServiceWorker() {
  return new Promise(() => null);
}