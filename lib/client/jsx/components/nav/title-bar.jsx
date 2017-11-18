import * as React from 'react'

export default class TitleBar extends React.Component{

  constructor(){

    super();
  }

  render(){
    
    return (

      <div id='title-menu'>

        <button className='title-menu-btn'>
            
          { 'metis' }
          <br />
          <span>

            { 'UPLOADER' }
          </span>
        </button>
        <img id='ucsf-logo' src='/img/ucsf_logo_dark.png' alt='' />
      </div>
    );
  }
}