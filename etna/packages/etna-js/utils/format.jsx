/*
 * Change a timestamp in to a human readable format.
 */

export const authorFormat = author => {
  const [ email, name ] = author.split('|');

  return { name, email }
}

export const dateFormat = (timestamp, missing='Unknown') => {
  if (timestamp == undefined || timestamp == null) return missing;

  let date = new Date(timestamp);
  let today = new Date();

  if (date.toDateString() == today.toDateString()) {
    let time = date.toLocaleString('en-us',{ hour: 'numeric', minute: 'numeric' });
    return `Today, ${time}`;
  } else {
    return date.toLocaleString('en-us', { month: 'short', day: 'numeric', year: 'numeric' });
  }
};

export const jsonFormat = obj => JSON.stringify(obj, null, 2);

/*
 * Change an integer of bytes into a human readable format.
 * http://stackoverflow.com/questions/10420352/converting-file-size-in-bytes-to-human-readable
 */

const PREFIXES = ['K','M','G','T','P','E','Z','Y']

export const byteFormat = (bytes, si = false, unit='B') => {
  let thresh = si ? 1000 : 1024;

  if (Math.abs(bytes) < thresh) return `${bytes} ${unit}`;

  let places = Math.min(
    PREFIXES.length,
    Math.floor(Math.log(bytes) / Math.log(thresh))
  );
  bytes = bytes / Math.pow(thresh, places);

  return `${bytes.toFixed(1)} ${PREFIXES[places - 1]}${si ? 'i' : ''}${unit}`;
};

export const snakeCase = (str) => {
  return str
    .split(/(?=[A-Z])/)
    .join('_')
    .toLowerCase();
};

export const camelCase = (str) => {
  return str
    .toLowerCase()
    .replace(/[^A-Za-z0-9]+([A-Za-z0-9])/g, (_, chr) => chr.toUpperCase());
};

export const capitalize = (str) => str.charAt(0).toUpperCase() + str.slice(1);
