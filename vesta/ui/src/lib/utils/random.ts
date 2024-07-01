export function getRandomItem<Item>(list: Item[]): Item {
    const randomIdx = Math.floor(Math.random() * list.length)
    return list[randomIdx]
}