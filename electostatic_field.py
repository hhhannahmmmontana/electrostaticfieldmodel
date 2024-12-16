import sys
from collections.abc import Callable
import tkinter as tk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
from tkinter import messagebox

canvas_margin = 1.25

def create_entry(frame: tk.Frame, var: tk.DoubleVar, row: int, column: int):
    tk.Entry(frame, textvariable=var).grid(row=row, column=column)

def create_label(frame: tk.Frame, text: str, row: int, column: int, fg_color: str = "white"):
    tk.Label(frame, text=text, fg=fg_color).grid(row=row, column=column, sticky=tk.W)

def create_input(frame: tk.Frame, text: str, var: tk.DoubleVar, row: int, column: int, text_color: str = "white"):
    create_label(frame, text, row, column, text_color)
    create_entry(frame, var, row, column + 1)

def create_button(frame: tk.Frame, text: str, func: Callable, row: int, column: int, padding_top: int, padding_bottom):
    tk.Button(
        frame,
        text=text,
        command=func
    ).grid(row=row, column=column, columnspan=2, pady=(padding_top, padding_bottom))

def create_listbox(frame: tk.Frame, row: int, column: int) -> tk.Listbox:
    listbox = tk.Listbox(frame, activestyle="none", width=30)
    listbox.grid(row=row, column=column, columnspan=2, sticky="nsew")
    return listbox

class ElectricFieldApp:
    def __init__(self, master: tk.Tk):
        self.root = master
        self.root.title("Моделирование электростатического поля")

        frame = tk.Frame(self.root, bg="black", padx=0)
        frame.pack(fill=tk.BOTH, expand=True)

        left_menu = tk.Frame(frame, bg="#333333", padx=10)
        left_menu.pack(fill=tk.Y, side=tk.LEFT)

        self.root.option_add("*Background", "#333333")
        self.root.option_add("*Foreground", "white")
        self.root.option_add("*Entry*Background", "#444444")
        self.root.option_add("*Entry*Relief", "flat")
        self.root.option_add("*Button*Background", "#444444")
        self.root.option_add("*Button*Relief", "flat")

        tk.Label(left_menu, text="Добавить заряд").grid(row=0, column=0, columnspan=2)

        self.x_var = tk.DoubleVar(value=0)
        create_input(left_menu, "x = ", self.x_var, 1, 0, "magenta")

        self.y_var = tk.DoubleVar(value=0)
        create_input(left_menu, "y = ", self.y_var, 2, 0, "magenta")

        self.q_var = tk.DoubleVar(value=0)
        create_input(left_menu, "q, Кл = ", self.q_var, 3, 0, "magenta")

        create_button(left_menu, "Добавить заряд", self.add_charge, 4, 0, 5, 10)

        self.dipole_x = tk.DoubleVar(value=0)
        create_input(left_menu, "x = ", self.dipole_x, 5, 0, "magenta")

        self.dipole_y = tk.DoubleVar(value=0)
        create_input(left_menu, "y = ", self.dipole_y, 6, 0, "magenta")

        self.dipole_length = tk.DoubleVar(value=0)
        create_input(left_menu, "l = ", self.dipole_length, 7, 0, "magenta")

        self.dipole_angle = tk.DoubleVar(value=0)
        create_input(left_menu, "\u03B8, \u00b0 = ", self.dipole_angle, 8, 0, "magenta")

        self.dipole_p = tk.DoubleVar(value=0)
        create_input(left_menu, "|p|, Кл = ", self.dipole_p, 9, 0, "magenta")

        create_button(left_menu, "Добавить диполь", self.add_dipole, 10, 0, 5, 10)

        create_button(left_menu, "Загрузить значения", self.visualize, 11, 0, 0, 10)
        self.charge_list = create_listbox(left_menu, 12, 0)
        create_button(left_menu, "Убрать значение", self.remove_selected, 13, 0, 10, 10)

        self.charges = []

        self.fig, self.ax = plt.subplots(figsize=(6, 5))
        self.canvas = FigureCanvasTkAgg(self.fig, master=frame)
        self.canvas_widget = self.canvas.get_tk_widget()

        self.create_canvas(frame)
        self.visualize()


    def create_canvas(self, frame: tk.Frame):
        self.canvas_widget.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        self.fig.patch.set_facecolor("#1A1A1A")
        self.ax.set_facecolor("#333333")

        self.ax.spines['top'].set_color('white')
        self.ax.spines['right'].set_color('white')
        self.ax.spines['left'].set_color('white')
        self.ax.spines['bottom'].set_color('white')

        self.ax.spines['top'].set_linewidth(2)
        self.ax.spines['right'].set_linewidth(2)
        self.ax.spines['left'].set_linewidth(2)
        self.ax.spines['bottom'].set_linewidth(2)

        self.ax.grid(True, color='white', linestyle='--', linewidth=0.5)

        self.ax.tick_params(axis='x', colors='white')
        self.ax.tick_params(axis='y', colors='white')

    def can_be_added(self, charge) -> bool:
        if charge['type'] == 'd' and charge['length'] == 0:
            messagebox.showerror("Не-а", "Длина диполя не может быть равна нулю!")
            return False
        # for c in self.charges:
        #     if charge['position'] == c['position']:
        #         messagebox.showerror("Не-а", "Диполь или заряд с такими координатами уже существует!")
        #         return False
        return True


    def add_dipole(self):
        x = self.dipole_x.get()
        y = self.dipole_y.get()
        l = self.dipole_length.get()
        theta = self.dipole_angle.get()
        p = self.dipole_p.get()

        charge = {"type": "d", "position": (x, y), "length": l, "angle": theta, "p": p}

        if not self.can_be_added(charge):
            return

        self.charges.append(charge)
        self.update_charge_list()

    def add_charge(self):
        x = self.x_var.get()
        y = self.y_var.get()
        q = self.q_var.get()

        charge = {"type": "c", "position": (x, y), "q": q}

        if not self.can_be_added(charge):
            return


        self.charges.append(charge)
        self.update_charge_list()

    def update_charge_list(self):
        self.charge_list.delete(0, tk.END)
        for charge in self.charges:
            text = (
                f"Заряд: q = {charge['q']}, x = {charge['position'][0]}, y = {charge['position'][1]}" if charge['type'] == 'c'
            else f"Диполь: l = {charge['length']}, \u03B8 = {charge['angle']}, |p| = {charge['p']}  x = {charge['position'][0]}, y = {charge['position'][1]}")
            self.charge_list.insert(tk.END, text)


    def remove_selected(self):
        selected = self.charge_list.curselection()
        if selected:
            index = selected[0]
            del self.charges[index]
            self.update_charge_list()

    def logic(self):
        min_x = sys.maxsize
        max_x = -(sys.maxsize - 1)
        min_y = sys.maxsize
        max_y = -(sys.maxsize - 1)
        for charge in self.charges:
            if charge['type'] == 'c':
                min_x = min(min_x, charge['position'][0] * canvas_margin)
                max_x = max(max_x, charge['position'][0] * canvas_margin)
                min_y = min(min_y, charge['position'][1] * canvas_margin)
                max_y = max(max_y, charge['position'][1] * canvas_margin)
            else:
                min_x = min(min_x,
                            charge['position'][0] + charge['length'] * np.cos(np.deg2rad(charge['angle'])),
                            charge['position'][0] - charge['length'] * np.cos(np.deg2rad(charge['angle'])))
                max_x = max(max_x,
                            charge['position'][0] + charge['length'] * np.cos(np.deg2rad(charge['angle'])),
                            charge['position'][0] - charge['length'] * np.cos(np.deg2rad(charge['angle'])))
                min_y = min(min_y,
                            charge['position'][1] + charge['length'] * np.sin(np.deg2rad(charge['angle'])),
                            charge['position'][1] - charge['length'] * np.sin(np.deg2rad(charge['angle'])))
                max_y = max(max_y,
                            charge['position'][1] + charge['length'] * np.sin(np.deg2rad(charge['angle'])),
                            charge['position'][1] - charge['length'] * np.sin(np.deg2rad(charge['angle'])))

        if max_x - min_x < 2:
            min_x -= 1
            max_x += 1
        if max_y - min_y < 2:
            min_y -= 1
            max_y += 1

        min_x = max_x - (max_x - min_x) * canvas_margin
        max_x = min_x + (max_x - min_x) * canvas_margin
        min_y = max_y - (max_y - min_y) * canvas_margin
        max_y = min_y + (max_y - min_y) * canvas_margin

        x = np.linspace(min_x, max_x, 800)
        y = np.linspace(min_y, max_y, 800)
        x_grid, y_grid = np.meshgrid(x, y)

        e_x = np.zeros(x_grid.shape)
        e_y = np.zeros(y_grid.shape)
        v = np.zeros(x_grid.shape)

        for charge in self.charges:
            if charge['type'] == 'c':
                q = charge["q"]
                cx, cy = charge["position"]
                dx = x_grid - cx
                dy = y_grid - cy
                r_squared = dx ** 2 + dy ** 2

                r_squared[r_squared < 1e-4] = 1e-4

                e_x += q * dx / r_squared ** (3 / 2)
                e_y += q * dy / r_squared ** (3 / 2)

                v += q / np.sqrt(r_squared)
            else:
                tx, ty = charge["position"]
                cx1 = tx + np.cos(np.deg2rad(charge["angle"])) * charge["length"]
                cy1 = ty + np.sin(np.deg2rad(charge["angle"])) * charge["length"]
                cx2 = tx - np.cos(np.deg2rad(charge["angle"])) * charge["length"]
                cy2 = ty - np.sin(np.deg2rad(charge["angle"])) * charge["length"]

                q1 = charge["p"] / charge["length"]
                q2 = -1 * charge["p"] / charge["length"]

                dx1 = x_grid - cx1
                dy1 = y_grid - cy1
                r1_squared = dx1 ** 2 + dy1 ** 2

                dx2 = x_grid - cx2
                dy2 = y_grid - cy2
                r2_squared = dx2 ** 2 + dy2 ** 2

                r1_squared[r1_squared < 1e-4] = 1e-4
                r2_squared[r2_squared < 1e-4] = 1e-4

                e_x += q1 * dx1 / r1_squared ** (3 / 2)
                e_y += q1 * dy1 / r1_squared ** (3 / 2)

                e_x += q2 * dx2 / r2_squared ** (3 / 2)
                e_y += q2 * dy2 / r2_squared ** (3 / 2)

                v += q1 / np.sqrt(r1_squared)
                v += q2 / np.sqrt(r2_squared)

        self.ax.streamplot(x_grid, y_grid, e_x, e_y, color="#67ff4d", linewidth=1, cmap='viridis', density=1.5)

        v_pos = v[v > 0]
        v_neg = v[v < 0]

        if len(v_pos) > 0:
            levels_pos = np.logspace(np.log10(np.nanmin(v_pos)), np.log10(np.nanmax(v_pos)), num=25)
            plt.contour(x, y, v, levels=levels_pos, colors="red", alpha=0.75)

        if len(v_neg) > 0:
            v_neg *= -1
            levels_neg = -1 * np.logspace(np.log10(np.nanmax(v_neg)), np.log10(np.nanmin(v_neg)), num=25)
            plt.contour(x, y, v, levels=levels_neg, colors="blue", alpha=0.75)

        for charge in self.charges:
            if (charge['type'] == 'c'):
                color = 'red' if charge["q"] > 0 else ('#0000FF' if (charge["q"] < 0) else 'white')
                self.ax.scatter(charge["position"][0], charge["position"][1], color=color, s=100, label=f"q = {charge['q']}", zorder=6)
            else:
                self.ax.plot(
                    [charge['position'][0] + np.cos(np.deg2rad(charge['angle'])) * charge['length'],
                    charge['position'][0] - np.cos(np.deg2rad(charge['angle'])) * charge['length']],
                    [charge['position'][1] + np.sin(np.deg2rad(charge['angle'])) * charge['length'],
                     charge['position'][1] - np.sin(np.deg2rad(charge['angle'])) * charge['length']],
                    linewidth=5,
                    color="magenta")

    def visualize(self):
        self.ax.clear()

        if len(self.charges) != 0:
            self.logic()

        self.ax.set_title("Моделирование электростатического поля", color="white")
        self.ax.set_xlabel("x", color="white")
        self.ax.set_ylabel("y", color="white")

        self.canvas.draw()
        self.canvas.get_tk_widget().update_idletasks()

def create_tkinter() -> tk.Tk:
    tmp = tk.Tk()
    tmp.geometry("800x600")
    tmp.option_add("*Font", "Arial 10")
    return tmp


if __name__ == "__main__":
    root = create_tkinter()
    app = ElectricFieldApp(root)
    root.mainloop()
