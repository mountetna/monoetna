import * as React from 'react'

export default class MenuBar extends React.Component{

  constructor(){

    super();
  }

  render(){
    
    return (

      <div id='nav-menu'>

        <div id='master-search-group'>

          <button id='master-search-button'>

            <span className='glyphicon glyphicon-search white-glyphicon'></span>
          </button>
        </div>
        
        {/*
        <button className='nav-menu-btn'>

          { 'ACTIVITY' }
        </button>
        <button className='nav-menu-btn'>

          { 'DOCS' }
        </button>
        <button className='nav-menu-btn'>

          { 'PLOT' }
        </button>
        <a className='nav-menu-btn' href='./experiments.html'>

          { 'EXPERIMENTS' }
        </a>
        */}

      </div>
    );
  } 
}