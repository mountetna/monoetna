// Mostly to use React.memo to try and
//   cache the images.
import React from 'react';

const ImageMemo = React.memo(({src, className, alt = ''}) => {
  return <img className={className} src={src} alt={alt} />;
});

export default ImageMemo;
