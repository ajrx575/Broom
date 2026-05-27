FROM python:3.10.19
WORKDIR /app
ADD broom.py .
ADD check_compression.py .
ADD process_gff.py .
ADD create_interm_ENS_files.py .
ADD create_interm_RS_files.py .
ADD create_interm_JGI_files.py .
ADD snpEff_exec.py .
ADD snpEff_v4_3t_core .
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
LABEL authors="Adia Redd"
WORKDIR /data
ENTRYPOINT ["python", "/app/broom.py"]
