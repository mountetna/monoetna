from helper import helper

@helper
def get_description():
    return str("""What is your favorite tool?
    
input: tool (str tool name)
output: manifesto (str of your matching manifesto)
""")

def func(tool = None):
    if tool == "sickle":
        return str('''A spectre is haunting Europe — the spectre of communism. All the powers of old Europe have entered into a holy alliance to exorcise this spectre: Pope and Tsar, Metternich and Guizot, French Radicals and German police-spies.\nWhere is the party in opposition that has not been decried as communistic by its opponents in power? Where is the opposition that has not hurled back the branding reproach of communism, against the more advanced opposition parties, as well as against its reactionary adversaries?

                # Two things result from this fact:

                # I. Communism is already acknowledged by all European powers to be itself a power.

                # II. It is high time that Communists should openly, in the face of the whole world, publish their views, their aims, their tendencies, and meet this nursery tale of the Spectre of Communism with a manifesto of the party itself.

                # To this end, Communists of various nationalities have assembled in London and sketched the following manifesto, to be published in the English, French, German, Italian, Flemish and Danish languages.''')
                
    elif tool == "musket":
        return str('''The unanimous Declaration of the thirteen united States of America, When in the Course of human events, it becomes necessary for one people to dissolve the political bands which have connected them with another, and to assume among the powers of the earth, the separate and equal station to which the Laws of Nature and of Nature's God entitle them, a decent respect to the opinions of mankind requires that they should declare the causes which impel them to the separation.''')
        
    elif tool == "bowling ball":
        return str('''The Dude abides''')
    
    else:
        return str('''I dont think nihlists have a manifesto''')