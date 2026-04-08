
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
    ys10 = ys[:,9]

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
    plt.figure()
    plt.plot(t, ys7, 'b')
    plt.xlim(0,days)
    plt.xlabel('Tempo (dias)')
    plt.ylabel('Concentração')
    plt.legend(["Cm"])
    plt.title('Células-tronco mesenquimais')
    plt.grid(True)
    plt.savefig(f'{path}cm.png')
    print(f'{path}cm.png')
    #plt.show()

    # Osteoblasts
    plt.figure()
    plt.plot(t, ys8, 'y')
    plt.legend(["Cb"])
    plt.xlim(0,days)
    plt.xlabel('Tempo (dias)')
    plt.ylabel('Concentração')
    plt.title('Osteoblastos')
    plt.grid(True)
    plt.savefig(f'{path}Cb.png')
    print(f'{path}Cb.png')
    
    # cartilage
    plt.figure()
    plt.plot(t, ys9, 'g', label="")
    #plt.plot(df["x"], df["y"], marker="o", linestyle="-", label="Modelo referência")
    plt.xlim(0,days)
    plt.xlabel('Time (days)')
    plt.ylabel('Concentration')
    plt.legend(["Mc"])
    plt.title('Cartilage')
    plt.grid(True)
    plt.legend()
    plt.savefig(f'{path}mc.png')
    print(f'{path}mc.png')
    
    
    #bone
    plt.figure()
    plt.plot(t, ys10, 'r', label="")
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
    
    
    plt.tight_layout()

def plot_estradiol(t, days, ys):
    ys11 = ys[:, 10]  # E2

    plt.figure()
    plt.plot(t, ys11, 'm', linewidth=2)
    plt.xlim(0, days)
    plt.xlabel('Tempo (dias)')
    plt.ylabel('Concentração')
    plt.title('Estradiol')
    plt.legend(['E2'])
    plt.grid(True)
    plt.savefig(f'{path}E2.png')
    print(f'{path}E2.png')
    plt.show()



