import React from 'react';
import ReactDOM from 'react-dom';
import { Provider } from 'react-redux';


function Vesta() {
    console.log("hello from vesta client lib!")
    // build the UI
    ReactDOM.render(
        <Provider>
            <></>
        </Provider>,
        document.getElementById('root')
    );
}

window.vesta = new Vesta();