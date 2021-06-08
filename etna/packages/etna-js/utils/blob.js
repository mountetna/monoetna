import downloadjs from 'downloadjs'

export const downloadBlob = ({ data, filename = 'file.dat', contentType = 'application/octet-stream' }) => {
  let blob = new Blob([data], {type: contentType});

  downloadjs(blob, filename, blob.type);
}

export const readTextFile = (accept) => (
  new Promise( (resolve, reject) => {
    let input = document.createElement("input");
    input.type = "file";
    input.accept = accept;
    input.addEventListener("change", e => {
      let file = e.target.files[0];


      let reader = new FileReader();
      reader.onload = () => resolve(reader.result);
      reader.onerror = () => reject(reader.error);
      reader.readAsText(file);
    });
    input.dispatchEvent(new MouseEvent("click")); 
  })
)
