# N-body simulator

Simulador gravitacional com modelos de vela solar por pressão de radiação e vela elétrica.

## Instalação

```bash
pip install -r requirements.txt
```

## Fluxo da e-sail

- `resources/electric_sail_probe.json`: exemplo de configuração da missão.
- `src/utils/schemas.py` e `src/model/electric_probe.py`: validam a entrada e criam a sonda.
- `src/services/n_body_service.py`: prepara os corpos e seleciona os dados do vento solar na data inicial.
- `src/n_body/n_body_solver.py`: integra a gravidade e soma a aceleração da vela.
- `src/solar_sail/electric_sail_dynamic.py`: calcula o empuxo e a aceleração da e-sail.

## Unidades da configuração

No JSON, `mass` é kg; `position` e `radius` são km; `velocity` é km/s. `radius` é o raio usado para colisões, não o raio dos fios.

Na e-sail, `N` é o número de fios, `L` é o comprimento de cada fio em km, `V` é a tensão em volts e `r_w` é o raio do fio em km. `phi_deg` e `theta_deg` são fornecidos em graus e convertidos para radianos ao criar a sonda.

O `mass` informado compõe o corpo usado pelo solver gravitacional. A aceleração propulsiva da e-sail usa a massa seca calculada pelo modelo (equações 75–76), não esse campo.

No fluxo normal, o serviço seleciona uma linha da tabela `src/atividade_solar.csv` na data inicial e mantém esses valores durante a integração. Chamadas diretas ao modelo, sem essa configuração, usam os valores padrão da classe.
