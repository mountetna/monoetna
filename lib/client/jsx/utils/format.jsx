/*
 * Change a timestamp in to a human readable format.
 */

export const userFormat = (author) => {
  let [ email, name ] = author.split(/\|/);
  return name;
}

export const dateFormat = (timestamp) => {
  if(timestamp == undefined || timestamp == null) return 'Unknown';

  let date = new Date(timestamp);
  let today = new Date();

  if (date.toDateString() == today.toDateString()) {
    let time = date.toLocaleString('en-us',{ hour: 'numeric', minute: 'numeric' });
    return `Today, ${time}`;
  } else {
    return date.toLocaleString('en-us', { month: 'short', day: 'numeric', year: 'numeric' });
  }
}

/*
 * Change an integer of bytes into a human readable format.  
 * http://stackoverflow.com/questions/10420352/converting-file-size-in-bytes-to-human-readable
 */

export const byteFormat = (bytes, si, bits = false) => {
  let thresh = si ? 1000 : 1024;
  if(Math.abs(bytes) < thresh){

      return bytes + ' B';
  }
  let units = si
      ? ['KB','MB','GB','TB','PB','EB','ZB','YB']
      : ['KiB','MiB','GiB','TiB','PiB','EiB','ZiB','YiB'];

  if (bits) units = ['Kbps','Mbps','Gbps'];

  let u = -1;
  do {
    bytes /= thresh;
    ++u;
  }
  while(Math.abs(bytes) >= thresh && u < units.length - 1);

  return bytes.toFixed(1)+' '+units[u];
}

export const snakeCase = (str) => {
  return str.split(/(?=[A-Z])/).join('_').toLowerCase();
}

export const camelCase = (str) => {
  return str.toLowerCase().replace(
    /[^A-Za-z0-9]+([A-Za-z0-9])/g, 
    (_, chr) =>  chr.toUpperCase()
  );
}
