
import matplotlib.pyplot as plt
import pandas as pd

path = "C:/Users/anacl/OneDrive/Área de Trabalho/Bone-estradiol-model-main/Bone-estradiol-model-main/output/"
#path = 'output/'

def plots(t,days,ys):
    
    ys1 = ys[:,0]
    ys2 = ys[:,1]
    ys3 = ys[:,2]
    ys4 = ys[:,3]
    ys5 = ys[:,4]
    ys6 = ys[:,5]
    ys7 = ys[:,6]
    ys8 = ys[:,7]
    ys9 = ys[:,8]
    #ys10 = ys[:,9]

    # Debris
    plt.figure()  
    plt.plot(t, ys1, 'r')
    plt.legend(["D"])
    plt.xlim(0,days)
    plt.xlabel('Tempo (dias)')
    plt.ylabel('Concentração')
    plt.title('Detritos')
    plt.legend(loc='upper right', fontsize=12)
    plt.grid(True)
    plt.savefig(f'{path}D.png')
    print(f'{path}D.png')
    #plt.show()


    # Macrofagos
    plt.figure()
    plt.plot(t, ys2, '--k', t,ys3, 'b', t,ys4, 'g')
    plt.legend(["M0", "M1", "M2"])
    plt.xlim(0,days)
    plt.xlabel('Tempo (dias)')
    plt.ylabel('Concentração')
    plt.title('Macrófagos')
    plt.grid(True)
    plt.savefig(f'{path}MoM1M2.png')
    print(f'{path}MoM1M2.png')
    #plt.show()

    # TNFalfa - pró-inflamatória
    plt.figure()
    plt.plot(t, ys5, 'g')
    plt.legend(["c1"])
    plt.xlim(0,days)
    plt.xlabel('Tempo (dias)')
    plt.ylabel('Concentração')
    plt.title('Citocinas pró-inflamatórias')
    plt.grid(True)
    plt.savefig(f'{path}c1.png')
    print(f'{path}c1.png')
    #plt.show()

    # IL10
    plt.figure()
    plt.plot(t, ys6, 'c')
    plt.legend(["c2"])
    plt.xlim(0,days)
    plt.xlabel('Tempo (dias)')
    plt.ylabel('Concentração')
    plt.title('Citocinas anti-inflamatórias')
    plt.grid(True)
    plt.savefig(f'{path}c2.png')
    print(f'{path}c2.png')
    
     # Mesenchymal stem cells
    #plt.figure()
    #plt.plot(t, ys7, 'b')
    #plt.xlim(0,days)
    #plt.xlabel('Tempo (dias)')
    #plt.ylabel('Concentração')
    #plt.legend(["Cm"])
    #plt.title('Células-tronco mesenquimais')
    #plt.grid(True)
    #plt.savefig(f'{path}cm.png')
    #print(f'{path}cm.png')
    #plt.show()

    # Osteoblasts
    plt.figure()
    plt.plot(t, ys7, 'y')
    plt.legend(["Cb"])
    plt.xlim(0,days)
    plt.xlabel('Tempo (dias)')
    plt.ylabel('Concentração')
    plt.title('Osteoblastos')
    plt.grid(True)
    plt.savefig(f'{path}Cb.png')
    print(f'{path}Cb.png')
    
 # Carregar o CSV usando delimitador ";" e decimal ","
    #df = pd.read_csv("article_mc.csv", delimiter=";", decimal=",")
    #dv = pd.read_csv("article_mb.csv", delimiter=";", decimal="," )
    
    # Renomear para clareza (opcional)
    #df.columns = ["x", "y"]
    #dv.columns = ["x", "y"]
    
    # cartilage
    #plt.figure()
    #plt.plot(t, ys9, 'g', label="")
    #plt.plot(df["x"], df["y"], marker="o", linestyle="-", label="Modelo referência")
    #plt.xlim(0,days)
    #plt.xlabel('Time (days)')
    #plt.ylabel('Concentration')
    #plt.legend(["Mc"])
    #plt.title('Cartilage')
    #plt.grid(True)
    #plt.legend()
    #plt.savefig(f'{path}mc.png')
    #print(f'{path}mc.png')
    
    
    #bone
    plt.figure()
    plt.plot(t, ys8, 'r', label="")
    #plt.plot(dv["x"], dv["y"], marker="o", linestyle="-", label="Modelo referência")
    plt.xlim(0,days)
    plt.xlabel('Time (days)')
    plt.ylabel('Concentration')
    plt.legend(["Mb"])
    plt.title('Bone formation')
    plt.grid(True)
    plt.legend()
    plt.savefig(f'{path}mb.png')
    print(f'{path}mb.png')
    
    
     #Estradiol
    plt.figure()
    plt.plot(t, ys9, 'r', label="")
    plt.xlim(0,days)
    plt.xlabel('Time (days)')
    plt.ylabel('Concentration')
    plt.legend(["Mb"])
    plt.title('Estradiol')
    plt.grid(True)
    plt.legend()
    plt.savefig(f'{path}E2.png')
    print(f'{path}E2.png')
    
    plt.tight_layout()

def plots_estradiol(t,days,ys, ys_min):
    # desacopla os resultados
    ys60 = ys[:,10] # estradiol E2 = 0.06
    ys60osso = ys[:,9]
    #ys60cart = ys[:,8]
    ys60c2 = ys[:,5]

    ys19 = ys_min[:,10] # estradiol E2 = 0.019
    ys19osso = ys_min[:,9]
    ys19cart = ys_min[:,8]
    ys19c2 = ys_min[:,5]

    d = {'tempo': t, 'osso': ys60osso}
    df = pd.DataFrame(data=d)
    # Write the DataFrame to a CSV file
    df.to_csv('osso_valormax.csv', index=False)
    print('Arquivo gerado') 
    #fazer para 0,19
    d2 = {'tempo': t, 'osso': ys19osso}
    df2 = pd.DataFrame(data=d2)
    print(df2)
    # Write the DataFrame to a CSV file
    df2.to_csv('osso_valormin.csv', index=False)
    
    fig, ax = plt.subplots(figsize=(5, 4), tight_layout=True)
    ax.plot(t, ys19, 'k', t, ys60, 'g')
    #ax.xlim(0,days)
    ax.legend(('E2 = 0.019', 'E2 = 0.060'), fontsize=14) #trocar
    ax.set_title('Estradiol', fontsize=14)
    ax.set_xlabel('Tempo (dias)', fontsize=14)
    ax.set_ylabel('Concentraçao', fontsize=14)
    ax.set_xlim(0,days)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f'{path}E2.png')

    #ys19osso = ys10

    fig, ax = plt.subplots(figsize=(5, 4), tight_layout=True)
    ax.plot(t, ys19osso, 'k', t, ys60osso, 'g')
    ax.legend(('E2 = 0.019', 'E2 = 0.060'), fontsize=14) #trocar
    ax.set_title('Osteocytes', fontsize=14)
    ax.set_xlabel('Time (days)', fontsize=14)
    ax.set_ylabel('Concentration (g/ml)', fontsize=14)
    ax.set_xlim(0,days)

    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f'{path}Osso_E2.png') 

    fig, ax = plt.subplots(figsize=(5, 4), tight_layout=True)
    ax.plot(t, ys19c2, 'k', t, ys60c2, 'g')
    ax.legend(('E2 = 0.019', 'E2 = 0.060'), fontsize=14) #trocar
    ax.set_title('IL-10', fontsize=14)
    ax.set_xlabel('Tempo (dias)', fontsize=14)
    ax.set_ylabel('Concentraçao', fontsize=14)
    ax.set_xlim(0,days)
    plt.grid(True)

    plt.tight_layout()
    plt.savefig(f'{path}IL-10_E2.png')

    #ys19osso = ys10

    #fig, ax = plt.subplots(figsize=(5, 4), tight_layout=True)
    #ax.plot(t, ys19cart, 'g', t, ys60cart, 'k')
    #ax.legend(('E2 = 0.019', 'E2 = 0.060'), fontsize=14) 
    #ax.set_title('Cartilagem', fontsize=14)
    #ax.set_xlabel('Tempo (dias)', fontsize=14)
    #ax.set_ylabel('Concentraçao', fontsize=14)
    #ax.set_xlim(0,days)
    #plt.grid(True)
    #plt.tight_layout()
    #plt.savefig(f'{path}Cartilagem_E2.png')



