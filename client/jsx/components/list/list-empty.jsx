import * as React from 'react'

export default class ListEmpty extends React.Component{

  constructor(){

    super();
  }

  render(){
    
    return (

      <div id='list-empty-group'>

        There are no files here.
      </div>
    );
  }
}