export function removeFromSet<T>(_set: Set<T>, items: T[]) {
    items.forEach(item => {
        _set.delete(item)
    })

    return _set
}