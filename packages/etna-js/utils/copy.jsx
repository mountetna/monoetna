export const copyText = (str) => {
  let text = document.createElement('textarea');
  text.value = str;
  text.style = {display: 'none'};

  document.body.appendChild(text);
  text.select();
  document.execCommand('copy');
  document.body.removeChild(text);
}
