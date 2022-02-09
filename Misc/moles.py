import time
from random import randint

def init():
  if randint(0,1) == 0:
    print("""
                               88             
                               88             
                               88             
88,dPYba,,adPYba,   ,adPPYba,  88  ,adPPYba, ,adPPYba,  
88P'   "88"    "8a a8"     "8a 88 a8P_____88 I8[    "" 
88      88      88 8b       d8 88 8PP\"\"\"\"\"\"\"  `'Y8ba,
88      88      88 "8a,   ,a8" 88 "8b,   ,aa aa    ]8I
88      88      88  `"YbbdP"'  88  `"Ybbd8"' `"YbbdP"'

MOLEcular dynamics Suite

                              _.-=-._                                    
                           o~`  '  > `.                                   
                           `.  ,       :                                   
                            `"-.__/    `.                                   
                               /       ::                                  
                              / .:    .:|                                 
                             :       .::!.                               
                            /'| ::  :::'!!                               
                          .:  :/' .::/  !!                               
                          :::/   :::'   !!                               
                          `:"::'''!!    !!                                 
                            /          :!!.                                
                           /     .-~-:  !!!                                
                          /:   :'        !!.                               
                         :::  :'          !!                               
                         |::  |        :!!!!                               
                         `::  :        !!!!'                               
                          |:. `:    .  '!!!                                
                          `::.  \   `::. !'                                
                           _.`::.\     ::                                  
                        .-~_____:~~    :'                                  
                        ~~~  .-'__..-~'           Loki                         
                             ~~~
""")
  else:
    print("""
                               88             
                               88             
                               88             
88,dPYba,,adPYba,   ,adPPYba,  88  ,adPPYba, ,adPPYba,  
88P'   "88"    "8a a8"     "8a 88 a8P_____88 I8[    "" 
88      88      88 8b       d8 88 8PP\"\"\"\"\"\"\"  `'Y8ba,
88      88      88 "8a,   ,a8" 88 "8b,   ,aa aa    ]8I
88      88      88  `"YbbdP"'  88  `"Ybbd8"' `"YbbdP"'

MOLEcular dynamics Suite

                              ,`""',
                              ;' ` ;
                              ;`,',;
                              ;' ` ;
                         ,,,  ;`,',;
                        ;,` ; ;' ` ;   ,',
                        ;`,'; ;`,',;  ;,' ;
                        ;',`; ;` ' ; ;`'`';
                        ;` '',''` `,',`',;
                         `''`'; ', ;`'`'
                              ;' `';
                              ;` ' ;
                              ;' `';
                              ;` ' ;
                              ; ',';
                              ;,' ';
                             \|/\|/;\|/           Pushi          
""")

    time.sleep(1.01)

if __name__=='__main__':
    init()