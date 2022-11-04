import { reducer }  from '../compose-identifier';

describe('reducer', () => {
  beforeEach(() => {
    store = mockStore({
      user: {
	email: 'janus@two-faces.org',
	name: 'Janus Bifrons',
        permissions: { 
          administration: {
            privileged: false,
            project_name: 'administration',
            role: 'editor'
          }
        }
      }
    });
  });

  it('sets some tokens', () => {
    const state = { tokens: [], seq: '' };
    const newState = reducer(
      state,
      {
        type: 'SET_TOKENS',
        tokens: [
          {
            "name": "PROJ",
            "label": "project",
            "values": {
              "LABORS": "The Twelve Labors of Hercules"
            }
          },
          {
            "name": "SEP",
            "label": "separator",
            "values": {
               "-": "# Separator"
             }
          },
          {
            "name": "LABOR",
            "label": "labor",
            "values": {
              "LION": "The Nemean Lion",
              "HYDRA": "The Lernean Hydra"
            }
          },
          {
            "name": "SEP",
            "label": "separator",
            "values": {
               "-": "# Separator"
             }
          },
          {
            "name": "VILL",
            "label": "village type",
            "values": {
              "V": "Village",
              "H": "Hamlet"
            }
          },
          {
                "name": "n",
                "label": "patient_counter"
          }
        ]
      }
    );
    expect( newState.height ).toEqual(4);
    expect( newState.regex ).toEqual("LABORS-(LION|HYDRA)-(V|H)[0-9]+");
    expect( newState.seq ).toEqual("LABORS-LABOR-VILLn");
    expect( newState.tokens.map(({assigned:_})=>_) ).toEqual([]);
    expect( newState.tokens.map(({from:_})=>_) ).toEqual([]);
    expect( newState.tokens.map(({to:_})=>_) ).toEqual([]);
    expect( newState.tokens.map(({filled:_})=>_) ).toEqual([]);
    newState.tokens.forEach(
      token => expect( Object.keys(token) ).toEqual( expect.arrayContaining(["assigned", "filled""from""height""label" "project", "name", "seq", "to", "type", "value", "values" ]) )
    )
  });
});
