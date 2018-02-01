/*
 * Change a timestamp in to a human readable format.
 */

export const dateFormat = (timestamp) => {
  timestamp = parseInt(timestamp);
  if(isNaN(timestamp) || timestamp < 0){
    timestamp = new Date().getTime() / 1000;
  }
  let dt = new Date(timestamp * 1000);
  let hr = dt.getHours() + 1;
  let mn = dt.getMinutes();
  let mo = dt.getMonth() + 1;
  let dy = dt.getDate();
  let yr = dt.getYear();

  return `${hr}:${mn}, ${mo}/${dy}`;
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
      ? ['kB','MB','GB','TB','PB','EB','ZB','YB']
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
