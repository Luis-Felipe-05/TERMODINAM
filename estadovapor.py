from tkinter import *
from tkinter import ttk
import bancopropriedadeinterface.rotinas_bancopropriedades as pd
import bancocpidealinterface.rotinas_banco_cp_ideal as cd
import modelos_termodinamicos.vapor as v
from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas
from reportlab.lib import colors
from pathlib import Path
#rotinas de cálculo
def exportar_pdf():
    # Pega o texto do resultado
    conteudo = resultado_box.get("1.0", END).strip()
    if not conteudo:
        return

    # Cria PDF
    #pdf_file = "dados.pdf"
    pdf_file = str(Path(__file__).parent.parent/"relatorios/dados_vapor.pdf")
    c = canvas.Canvas(pdf_file, pagesize=A4)
    largura, altura = A4

    # Fundo escuro
    c.setFillColorRGB(0.04, 0.06, 0.12)  # equivalente a #0a0f1f
    c.rect(0, 0, largura, altura, stroke=0, fill=1)

    # Título
    c.setFont("Helvetica-Bold", 20)
    c.setFillColorRGB(0.0, 0.96, 0.83)  # verde neon #00f5d4
    c.drawString(50, altura - 80, "RELATÓRIO")

    # Linha divisória neon
    c.setStrokeColorRGB(0.62, 0.30, 0.86)  # roxo neon #9d4edd
    c.setLineWidth(2)
    c.line(50, altura - 90, largura - 50, altura - 90)

    # Conteúdo em estilo console futurista
    y = altura - 120
    c.setFont("Courier", 11)
    c.setFillColor(colors.whitesmoke)

    for linha in conteudo.split("\n"):
        c.drawString(60, y, linha)
        y -= 15
        if y < 50:  # quebra de página
            c.showPage()
            c.setFillColorRGB(0.04, 0.06, 0.12)
            c.rect(0, 0, largura, altura, stroke=0, fill=1)
            c.setFont("Courier", 11)
            c.setFillColor(colors.whitesmoke)
            y = altura - 50

    c.save()
def delete():
    resultado_box.delete(1.0,END)
def resultado():
    if cb_eos.get() == "PR":
        t = float(temperatra_entry.get())
        p = float(pressao_entry.get())
        for tc in pd.consultar("SELECT R_TC FROM propriedades WHERE T_COMPOSTO = '"+cb_composto.get()+"'")[0]:
                tc
        for pc in pd.consultar("SELECT R_PC FROM propriedades WHERE T_COMPOSTO = '"+cb_composto.get()+"'")[0]:
                pc
        for w in pd.consultar("SELECT R_W FROM propriedades WHERE T_COMPOSTO = '"+cb_composto.get()+"'")[0]:
                w
        for a in cd.consultar("SELECT R_A FROM cp_gas_ideal WHERE T_COMPOSTO = '"+cb_composto.get()+"'")[0]:
                a
        for b in cd.consultar("SELECT R_B FROM cp_gas_ideal WHERE T_COMPOSTO = '"+cb_composto.get()+"'")[0]:
                b
        for c in cd.consultar("SELECT R_C FROM cp_gas_ideal WHERE T_COMPOSTO = '"+cb_composto.get()+"'")[0]:
                c
        for d in cd.consultar("SELECT R_D FROM cp_gas_ideal WHERE T_COMPOSTO = '"+cb_composto.get()+"'")[0]:
                d
        tol = 0.0001
        itmax = 100
        componente_class = v.Prop(p,t,float(pc),float(tc),float(w),itmax,tol)
        h = componente_class.entalpia_absolute_pr(float(a),float(b),float(c),float(d))
        s = componente_class.entropia_absolute_pr(float(a),float(b),float(c),float(d))
        ah = componente_class.helmholtz_absolute_pr(float(a),float(b),float(c),float(d))
        g = componente_class.gibbs_absolute_pr(float(a),float(b),float(c),float(d))
        f = componente_class.coeficiente_fug_pr()[1]
        cf = componente_class.coeficiente_fug_pr()[0]
        msg = f"ENTALPIA(KJ/KMOL): {h}\nENTROPIA(KJ/KMOL.K): {s}\nHELMOTZ(KJ/KMOL): {ah}\nGIBBS(KJ/KMOL): {g}\nFUGACIDADE(BAR): {f}\nCOEFICIENTE_FUGACIDADE: {cf}\nTEMPERATURA(K): {t}  PRESSÃO(bar):  {p}  MODELO: {cb_eos.get()} COMPONENTE: {cb_composto.get()}\nAUTOR: LUÍS FELIPE VILLALBA PEREIRA"
        #print(msg)
        resultado_box.insert('end',msg)
    if cb_eos.get() == "SRK":
        t = float(temperatra_entry.get())
        p = float(pressao_entry.get())
        for tc in pd.consultar("SELECT R_TC FROM propriedades WHERE T_COMPOSTO = '"+cb_composto.get()+"'")[0]:
                tc
        for pc in pd.consultar("SELECT R_PC FROM propriedades WHERE T_COMPOSTO = '"+cb_composto.get()+"'")[0]:
                pc
        for w in pd.consultar("SELECT R_W FROM propriedades WHERE T_COMPOSTO = '"+cb_composto.get()+"'")[0]:
                w
        for a in cd.consultar("SELECT R_A FROM cp_gas_ideal WHERE T_COMPOSTO = '"+cb_composto.get()+"'")[0]:
                a
        for b in cd.consultar("SELECT R_B FROM cp_gas_ideal WHERE T_COMPOSTO = '"+cb_composto.get()+"'")[0]:
                b
        for c in cd.consultar("SELECT R_C FROM cp_gas_ideal WHERE T_COMPOSTO = '"+cb_composto.get()+"'")[0]:
                c
        for d in cd.consultar("SELECT R_D FROM cp_gas_ideal WHERE T_COMPOSTO = '"+cb_composto.get()+"'")[0]:
                d
        tol = 0.0001
        itmax = 100
        componente_class = v.Prop(p,t,float(pc),float(tc),float(w),itmax,tol)
        h = componente_class.entalpia_absolute_srk(float(a),float(b),float(c),float(d))
        s = componente_class.entropia_absolute_srk(float(a),float(b),float(c),float(d))
        ah = componente_class.helmholtz_absolute_srk(float(a),float(b),float(c),float(d))
        g = componente_class.gibbs_absolute_srk(float(a),float(b),float(c),float(d))
        f = componente_class.coeficiente_fug_srk()[1]
        cf = componente_class.coeficiente_fug_srk()[0]
        msg = f"ENTALPIA(KJ/KMOL): {h}\nENTROPIA(KJ/KMOL.K): {s}\nHELMOTZ(KJ/KMOL): {ah}\nGIBBS(KJ/KMOL): {g}\nFUGACIDADE(BAR): {f}\nCOEFICIENTE_FUGACIDADE: {cf}\nTEMPERATURA(K): {t}  PRESSÃO(bar):  {p}  MODELO: {cb_eos.get()} COMPONENTE: {cb_composto.get()}\nAUTOR: LUÍS FELIPE VILLALBA PEREIRA"
        #print(msg)
        resultado_box.insert('end',msg)
    if cb_eos.get() == "RK":
        t = float(temperatra_entry.get())
        p = float(pressao_entry.get())
        for tc in pd.consultar("SELECT R_TC FROM propriedades WHERE T_COMPOSTO = '"+cb_composto.get()+"'")[0]:
                tc
        for pc in pd.consultar("SELECT R_PC FROM propriedades WHERE T_COMPOSTO = '"+cb_composto.get()+"'")[0]:
                pc
        for w in pd.consultar("SELECT R_W FROM propriedades WHERE T_COMPOSTO = '"+cb_composto.get()+"'")[0]:
                w
        for a in cd.consultar("SELECT R_A FROM cp_gas_ideal WHERE T_COMPOSTO = '"+cb_composto.get()+"'")[0]:
                a
        for b in cd.consultar("SELECT R_B FROM cp_gas_ideal WHERE T_COMPOSTO = '"+cb_composto.get()+"'")[0]:
                b
        for c in cd.consultar("SELECT R_C FROM cp_gas_ideal WHERE T_COMPOSTO = '"+cb_composto.get()+"'")[0]:
                c
        for d in cd.consultar("SELECT R_D FROM cp_gas_ideal WHERE T_COMPOSTO = '"+cb_composto.get()+"'")[0]:
                d
        tol = 0.0001
        itmax = 100
        componente_class = v.Prop(p,t,float(pc),float(tc),float(w),itmax,tol)
        h = componente_class.entalpia_absolute_rk(float(a),float(b),float(c),float(d))
        s = componente_class.entropia_absolute_rk(float(a),float(b),float(c),float(d))
        ah = componente_class.helmholtz_absolute_rk(float(a),float(b),float(c),float(d))
        g = componente_class.gibbs_absolute_rk(float(a),float(b),float(c),float(d))
        f = componente_class.coeficiente_fug_rk()[1]
        cf = componente_class.coeficiente_fug_rk()[0]
        msg = f"ENTALPIA(KJ/KMOL): {h}\nENTROPIA(KJ/KMOL.K): {s}\nHELMOTZ(KJ/KMOL): {ah}\nGIBBS(KJ/KMOL): {g}\nFUGACIDADE(BAR): {f}\nCOEFICIENTE_FUGACIDADE: {cf}\nTEMPERATURA(K): {t}  PRESSÃO(bar):  {p}  MODELO: {cb_eos.get()} COMPONENTE: {cb_composto.get()}\nAUTOR: LUÍS FELIPE VILLALBA PEREIRA"
        #print(msg)
        resultado_box.insert('end',msg)
    if cb_eos.get() == "VDW":
        t = float(temperatra_entry.get())
        p = float(pressao_entry.get())
        for tc in pd.consultar("SELECT R_TC FROM propriedades WHERE T_COMPOSTO = '"+cb_composto.get()+"'")[0]:
                tc
        for pc in pd.consultar("SELECT R_PC FROM propriedades WHERE T_COMPOSTO = '"+cb_composto.get()+"'")[0]:
                pc
        for w in pd.consultar("SELECT R_W FROM propriedades WHERE T_COMPOSTO = '"+cb_composto.get()+"'")[0]:
                w
        for a in cd.consultar("SELECT R_A FROM cp_gas_ideal WHERE T_COMPOSTO = '"+cb_composto.get()+"'")[0]:
                a
        for b in cd.consultar("SELECT R_B FROM cp_gas_ideal WHERE T_COMPOSTO = '"+cb_composto.get()+"'")[0]:
                b
        for c in cd.consultar("SELECT R_C FROM cp_gas_ideal WHERE T_COMPOSTO = '"+cb_composto.get()+"'")[0]:
                c
        for d in cd.consultar("SELECT R_D FROM cp_gas_ideal WHERE T_COMPOSTO = '"+cb_composto.get()+"'")[0]:
                d
        tol = 0.0001
        itmax = 100
        componente_class = v.Prop(p,t,float(pc),float(tc),float(w),itmax,tol)
        h = componente_class.entalpia_absolute_vw(float(a),float(b),float(c),float(d))
        s = componente_class.entropia_absolute_vw(float(a),float(b),float(c),float(d))
        ah = componente_class.helmholtz_absolute_vw(float(a),float(b),float(c),float(d))
        g = componente_class.gibbs_absolute_vw(float(a),float(b),float(c),float(d))
        f = componente_class.coeficiente_fug_vw()[1]
        cf = componente_class.coeficiente_fug_vw()[0]
        msg = f"ENTALPIA(KJ/KMOL): {h}\nENTROPIA(KJ/KMOL.K): {s}\nHELMOTZ(KJ/KMOL): {ah}\nGIBBS(KJ/KMOL): {g}\nFUGACIDADE(BAR): {f}\nCOEFICIENTE_FUGACIDADE: {cf}\nTEMPERATURA(K): {t}  PRESSÃO(bar):  {p}  MODELO: {cb_eos.get()} COMPONENTE: {cb_composto.get()}\nAUTOR: LUÍS FELIPE VILLALBA PEREIRA"
        #print(msg)
        resultado_box.insert('end',msg)

#--------------------------------------------------------------------------
root = Tk()
root.title("VAPOR/COMPOSTOS PUROS")
root.geometry("600x450")
root.configure(bg="#0a0f1f")  # fundo futurista escuro

# Paleta de cores
bg_color = "#0a0f1f"
label_color = "#00f5d4"   # verde neon
accent_color = "#9d4edd"  # roxo neon
entry_bg = "#111827"
entry_fg = "#e0e0e0"

# Deixa a janela responsiva
for i in range(3):
    root.columnconfigure(i, weight=1)
for i in range(9):
    root.rowconfigure(i, weight=1)

# Estilo para Combobox
root.option_add("*TCombobox*Listbox*Background", entry_bg)
root.option_add("*TCombobox*Listbox*Foreground", entry_fg)
root.option_add("*TCombobox*Listbox*selectBackground", accent_color)
root.option_add("*TCombobox*Listbox*selectForeground", "white")

# ---------- Widgets ----------
frameprincipal = Label(root, text="PROPRIEDADES",
                       bg=bg_color, fg=label_color,
                       font=("Segoe UI", 12, "bold"))
frameprincipal.grid(row=0, column=0, columnspan=3, pady=5)

selec_composto = Label(root, text="SELECIONE O COMPOSTO",
                       bg=bg_color, fg=label_color,
                       font=("Segoe UI", 10))
selec_composto.grid(row=1, column=0, sticky="w", padx=10)

query = "SELECT T_COMPOSTO FROM propriedades"
compostos = pd.consultar(query)

eq_estado = ['PR', 'SRK', 'RK', 'VDW']

cb_composto = ttk.Combobox(root, values=compostos, font=("Segoe UI", 10))
cb_composto.grid(row=1, column=1, columnspan=2, sticky="we", padx=10, pady=5)

temperatura = Label(root, text="TEMPERATURA (K):",
                    bg=bg_color, fg=label_color,
                    font=("Segoe UI", 10))
temperatura.grid(row=2, column=0, sticky="w", padx=10)

temperatra_entry = Entry(root, bg=entry_bg, fg=entry_fg,
                         insertbackground="white", relief="flat",
                         font=("Consolas", 10))
temperatra_entry.grid(row=2, column=1, sticky="we", padx=10, pady=5)

pressao = Label(root, text="PRESSÃO (bar):",
                bg=bg_color, fg=label_color,
                font=("Segoe UI", 10))
pressao.grid(row=3, column=0, sticky="w", padx=10)

pressao_entry = Entry(root, bg=entry_bg, fg=entry_fg,
                      insertbackground="white", relief="flat",
                      font=("Consolas", 10))
pressao_entry.grid(row=3, column=1, sticky="we", padx=10, pady=5)

eos = Label(root, text="SELECIONE O MODELO TERMODINÂMICO",
            bg=bg_color, fg=label_color,
            font=("Segoe UI", 10))
eos.grid(row=4, column=0, sticky="w", padx=10)

cb_eos = ttk.Combobox(root, values=eq_estado, font=("Segoe UI", 10))
cb_eos.grid(row=4, column=1, columnspan=2, sticky="we", padx=10, pady=5)

# ---------- Botões ----------
btn_calcular = Button(root, text="CALCULAR",
                      bg=accent_color, fg="white",
                      font=("Segoe UI", 10, "bold"),
                      activebackground=label_color, activeforeground="black",
                      relief="flat",
                      command=resultado)
btn_calcular.grid(row=5, column=0, padx=10, pady=10, sticky="we")

btn_limpar = Button(root, text="LIMPAR",
                    bg=label_color, fg="black",
                    font=("Segoe UI", 10, "bold"),
                    activebackground=accent_color, activeforeground="white",
                    relief="flat",
                    command=delete)
btn_limpar.grid(row=5, column=1, padx=10, pady=10, sticky="we")

btn_exportar = Button(root, text="EXPORTAR PDF",
                      bg="#ff0054", fg="white",
                      font=("Segoe UI", 10, "bold"),
                      activebackground=label_color, activeforeground="black",
                      relief="flat",
                      command=exportar_pdf)
btn_exportar.grid(row=5, column=2, padx=10, pady=10, sticky="we")

# ---------- Campo de resultado ----------
resultado_label = Label(root, text="RESULTADO:",
                        bg=bg_color, fg=label_color,
                        font=("Segoe UI", 10, "bold"))
resultado_label.grid(row=6, column=0, sticky="w", padx=10)

resultado_box = Text(root, height=6,
                     bg=entry_bg, fg=label_color,  # verde neon no texto
                     insertbackground="white", relief="flat",
                     font=("Consolas", 10))
resultado_box.grid(row=7, column=0, columnspan=3, sticky="nsew", padx=10, pady=5)

# ---------- Iniciar ----------
root.mainloop()

#resultado()