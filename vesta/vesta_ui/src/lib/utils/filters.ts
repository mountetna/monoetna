import { DATA_TYPES } from '@/lib/fixtures';

function random(seed) {
  let x = Math.sin(seed) * 10000;
  return x - Math.floor(x);
}
const dummyProjectCount = (project, dataType) => {
  // random number between 0 and 50
  let key = project + dataType;
  let keynum = 0;
  for (let i = 0; i < key.length; i++) {
    keynum += key.charCodeAt(i);
  }
  return parseInt(random(keynum) * 250);
}

export const projectDataTypes = project => [ ...new Set(
  project.dataTypes.map(
    modelName => Object.keys(DATA_TYPES).find(
      dt => DATA_TYPES[dt].includes(modelName)
    )
  ).filter(_=>_)
) ]
